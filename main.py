from flask import Flask, render_template, request, redirect, url_for, flash, send_file, session, abort
import os
from flask_cors import CORS
from tempfile import NamedTemporaryFile

from datetime import timedelta

from algorithms.Protein3D import StructureVisualisation

app = Flask(__name__)
CORS(app)
app.secret_key = 'your_secret_key'
app.config['UPLOAD_FOLDER'] = 'pdb'
app.config['ALLOWED_EXTENSIONS'] = {'pdb', 'pdb1'}

app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=1)


def allowed_file(filename):
    _, file_extension = os.path.splitext(filename)
    return file_extension in app.config['ALLOWED_EXTENSIONS']


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/pdb/structure.pdb")
def serve_pdb_file():
    return send_file("pdb/structure.pdb", as_attachment=True)


@app.route("/upload_structure", methods=['POST'])
def upload_structure():

    if 'file' not in request.files:
        abort(400)

    file = request.files['file']

    if not file:
        abort(400)

    try:
        sid = id(session)
        temp_file = NamedTemporaryFile(delete=False, mode='w+')
        temp_file.write(file.read().decode('utf-8'))
        temp_file.seek(0)
        structure = StructureVisualisation(sid, temp_file)
        session.structure = structure
        temp_file.close()
    except Exception:
        abort(400)

    return "200"



if __name__ == "__main__":
    app.run(host="localhost", port=8000, debug=True)
