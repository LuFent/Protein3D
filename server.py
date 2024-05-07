from flask import Flask, render_template, request, send_file, session, abort, jsonify, send_from_directory
import os
import atexit
from flask_cors import CORS
from tempfile import NamedTemporaryFile
import secrets
from algorithms.Algorithm import Algorithm, AlgorithmsStorage
from datetime import timedelta

from algorithms.StructureVisualisation import StructureVisualisation, NotDisordered


server = Flask(__name__)

STRUCTURES_STORE = {}
TEMP_FILES = set()


def requires_structure(func):
    def wrserverer(*args, **kwargs):
        if 'sid' not in session:
            return jsonify({'error': 'Session does not have the structure'}), 400
        return func(*args, **kwargs)

    wrserverer.__name__ = func.__name__
    return wrserverer


def allowed_file(filename):
    _, file_extension = os.path.splitext(filename)
    return file_extension in server.config['ALLOWED_EXTENSIONS']


def clean_tempfiles_dir():
    directory = server.config['TEMPFILE_DIR']
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)


@server.route("/")
def index():
    return render_template("index.html")


@server.route("/app_presets", methods=['GET'])
def get_app_presets():
    return {
        "algorithms": server.config["ALGORITHMS"].get_all_algorithms(),
    }


@server.route("/delete_structure_file", methods=['DELETE'])
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

        os.remove(os.path.join(server.config["TEMPFILE_DIR"], file_name))
    except Exception:
        pass
    return "ok", 200


@server.route("/upload_structure", methods=['POST'])
def upload_structure():
    if 'file' not in request.files:
        abort(400)

    file = request.files['file']
    filename = file.filename
    if not file:
        abort(400)

    if session.get("sid", False):
        session_id = session["sid"]
        if STRUCTURES_STORE.get(session_id, False):
            del STRUCTURES_STORE[session_id]
    else:
        session_id = secrets.token_hex(16)
        session['sid'] = session_id

    temp_file = NamedTemporaryFile(delete=True, mode='w+', dir=server.config['TEMPFILE_DIR'])
    temp_file.write(file.read().decode('utf-8'))
    temp_file.seek(0)
    structure = StructureVisualisation.get_structure(session_id, temp_file)

    temp_file.close()
    cleaned_temp_file = NamedTemporaryFile(delete=False, mode='w+', dir=server.config['TEMPFILE_DIR'])
    StructureVisualisation.clean_file(structure, cleaned_temp_file)
    cleaned_temp_file.seek(0)
    cleaned_structure = StructureVisualisation(session_id,
                                               cleaned_temp_file,
                                               server.config["ALGORITHMS"].copy(),
                                               filename)

    STRUCTURES_STORE[session_id] = cleaned_structure
    response = send_file(cleaned_temp_file.name, as_attachment=True)
    response.headers["Content-Disposition"] = os.path.basename(cleaned_temp_file.name)
    cleaned_temp_file.close()
    TEMP_FILES.add(cleaned_temp_file.name)

    return response


@server.route("/exec_algorithm", methods=['POST'])
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

    try:
        mask = structure.execute_algorithm(alg).serialize()
    except Exception:
        return jsonify({'error': 'Something went wrong'}), 400

    return jsonify(mask), 200


def run_server(*algs):
    algorithm_storage = AlgorithmsStorage()
    for alg in algs:
        algorithm_storage.add_algorithm(alg())
    server.config["ALGORITHMS"] = algorithm_storage
    server.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=True)
    atexit.register(clean_tempfiles_dir)
    CORS(server, supports_credentials=True, expose_headers='Content-Disposition')
    server.secret_key = 'your_secret_key'
    server.config['UPLOAD_FOLDER'] = 'pdb'
    server.config['ALLOWED_EXTENSIONS'] = {'pdb', 'pdb1'}
    server.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=1)
    server.config['TEMPFILE_DIR'] = "tempfiles"
    os.makedirs(server.config['TEMPFILE_DIR'], exist_ok=True)
    server.run(host="localhost", port=8000)
