from flask import Flask, render_template, request, redirect, url_for, flash, send_file, session, abort, jsonify
import os
import sys
import atexit
from flask_cors import CORS
from tempfile import NamedTemporaryFile
import secrets
import shutil
from Bio.PDB import Selection, PDBIO, Select
from datetime import timedelta

from algorithms.Protein3D import StructureVisualisation, NotDisordered


app = Flask(__name__)
CORS(app, supports_credentials=True, expose_headers='Content-Disposition')
app.secret_key = 'your_secret_key'
app.config['UPLOAD_FOLDER'] = 'pdb'
app.config['ALLOWED_EXTENSIONS'] = {'pdb', 'pdb1'}
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=1)
app.config['TEMPFILE_DIR'] = "tempfiles"
app.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=True)
STRUCTURES_STORE = {}
TEMP_FILES = set()


os.makedirs(app.config['TEMPFILE_DIR'], exist_ok=True)


def requires_structure(func):
    def wrapper(*args, **kwargs):
        if 'sid' not in session:
            return jsonify({'error': 'Session does not have the structure'}), 400
        return func(*args, **kwargs)

    wrapper.__name__ = func.__name__
    return wrapper


def allowed_file(filename):
    _, file_extension = os.path.splitext(filename)
    return file_extension in app.config['ALLOWED_EXTENSIONS']


def clean_tempfiles_dir():
    directory = app.config['TEMPFILE_DIR']
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)


@app.route("/")
def index():
    return render_template("index.html")



@app.route("/delete_structure_file", methods=['DELETE'])
def delete_structure_file():
    if request.is_json:
        json_data = request.json
        if 'file_name' in json_data:
            file_name = json_data['file_name']
        else:
            return jsonify({'error': 'JSON data does not contain the "file_name" key'}), 400
    else:
        return jsonify({'error': 'Request data is not JSON'}), 400

    try:

        os.remove(os.path.join(app.config["TEMPFILE_DIR"], file_name))
    except Exception:
        pass
    return "ok", 200


@app.route("/upload_structure", methods=['POST'])
def upload_structure():
    if 'file' not in request.files:
        abort(400)

    file = request.files['file']

    if not file:
        abort(400)

    if session.get("sid", False):
        session_id = session["sid"]
        if STRUCTURES_STORE.get(session_id, False):
            del STRUCTURES_STORE[session_id]
    else:
        session_id = secrets.token_hex(16)
        session['sid'] = session_id

    temp_file = NamedTemporaryFile(delete=True, mode='w+', dir=app.config['TEMPFILE_DIR'])
    temp_file.write(file.read().decode('utf-8'))
    temp_file.seek(0)
    structure = StructureVisualisation.get_structure(session_id, temp_file)
    temp_file.close()
    cleaned_temp_file = NamedTemporaryFile(delete=False, mode='w+', dir=app.config['TEMPFILE_DIR'])
    StructureVisualisation.clean_file(structure, cleaned_temp_file)
    cleaned_temp_file.seek(0)
    cleaned_structure = StructureVisualisation(session_id, cleaned_temp_file)
    STRUCTURES_STORE[session_id] = cleaned_structure
    response = send_file(cleaned_temp_file.name, as_attachment=True)
    response.headers["Content-Disposition"] = os.path.basename(cleaned_temp_file.name)
    cleaned_temp_file.close()
    TEMP_FILES.add(cleaned_temp_file.name)

    return response


@app.route("/exec_algorithm", methods=['POST'])
@requires_structure
def exec_algorithm():
    if request.is_json:
        json_data = request.json
        if 'alg' in json_data:
            alg = json_data['alg']
        else:
            return jsonify({'error': 'JSON data does not contain the "alg" key'}), 400
    else:
        return jsonify({'error': 'Request data is not JSON'}), 400

    structure = STRUCTURES_STORE[session["sid"]]
    mask = structure.execute_algorithm(alg)
    if mask is None:
        return jsonify({'error': f'No such algorithm {alg}'}), 400
    return jsonify(mask), 200



atexit.register(clean_tempfiles_dir)

if __name__ == "__main__":
    app.run(host="localhost", port=8000, debug=True)
