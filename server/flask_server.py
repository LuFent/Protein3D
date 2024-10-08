from pprint import pprint

from flask import Flask, render_template, request, send_file, abort, jsonify, send_from_directory
import os
import atexit
from flask_cors import CORS
from tempfile import NamedTemporaryFile
import secrets
from server.Algorithm.Algorithm import Algorithm, AlgorithmsStorage
from datetime import timedelta
from server.Algorithm.StructureVisualisation import StructureVisualisation, NotDisordered
from uuid import uuid4
from pprint import pprint
from functools import wraps


STRUCTURES_STORE = {}

server = Flask(__name__)
PORT = 8000

TEMP_FILES = set()



def create_session() -> str:
    session_id = str(uuid4())
    STRUCTURES_STORE[session_id] = None
    return session_id


def get_session_data(session_id: str) -> dict:
    return STRUCTURES_STORE.get(session_id)


def save_session_data(session_id: str, data: dict):
    STRUCTURES_STORE[session_id] = data


def get_structure_visualisation(file, session_id, filename):
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

    return cleaned_structure, cleaned_temp_file


def require_session_id(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        session_id = request.args.get('session_id')
        if not session_id:
            return jsonify({"error": "Missing session_id"}), 400
        if session_id not in STRUCTURES_STORE:
            return jsonify({"error": "No such session_id"}), 404
        return f(*args, **kwargs)
    return decorated_function


def require_structure(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        session_id = request.args.get('session_id')
        if not session_id:
            return jsonify({"error": "Missing session_id"}), 400
        if not STRUCTURES_STORE.get(session_id, None):
            return jsonify({"error": "No such session_id"}), 404
        return f(*args, **kwargs)
    return decorated_function



def allowed_file(filename):
    _, file_extension = os.path.splitext(filename)
    return file_extension in server.config['ALLOWED_EXTENSIONS']


def clean_tempfiles_dir():
    directory = server.config['TEMPFILE_DIR']
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)


#Корневой адрес
@server.route("/")
def index():
    return render_template("index.html")


#Все доступные алгоритмы
@server.route("/app_presets", methods=['GET'])
def get_app_presets():
    session_id = create_session()
    return {
        "algorithms": server.config["ALGORITHMS"].get_all_algorithms(),
        "session_id": session_id
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


#Загрузка структуры
@server.route("/upload_structure", methods=['POST'])
@require_session_id
def upload_structure():
    session_id = request.args.get('session_id')
    if 'file' not in request.files:
        abort(400)

    file = request.files['file']
    filename = file.filename
    if not file:
        abort(400)

    if STRUCTURES_STORE.get(session_id, False):
        del STRUCTURES_STORE[session_id]

    structure_visualisation, cleaned_temp_file = get_structure_visualisation(file, session_id, filename)
    STRUCTURES_STORE[session_id] = structure_visualisation
    response = send_file(cleaned_temp_file.name, as_attachment=True)
    cleaned_temp_file.close()
    TEMP_FILES.add(cleaned_temp_file.name)

    return response


#Выполнить алгоритм
@server.route("/exec_algorithm", methods=['POST'])
@require_structure
def exec_algorithm():
    session_id = request.args.get('session_id')
    if request.is_json:
        json_data = request.json
        if 'alg' in json_data:
            alg = json_data['alg']
        else:
            return jsonify({'error': 'JSON data does not contain the "alg" key'}), 400
    else:
        return jsonify({'error': 'Request data is not JSON'}), 400

    structure = STRUCTURES_STORE[session_id]

    try:
        mask = structure.execute_algorithm(alg).serialize()
    except Exception:
        return jsonify({'error': 'Something went wrong'}), 400
    return jsonify(mask), 200


def run_server(*algs, port=8000):
    algorithm_storage = AlgorithmsStorage()
    for alg in algs:
        algorithm_storage.add_algorithm(alg())
    server.config["ALGORITHMS"] = algorithm_storage
    server.config['SESSION_COOKIE_SAMESITE'] = 'None'
    server.config['SESSION_COOKIE_SECURE'] = False
    server.config['SESSION_COOKIE_HTTPONLY'] = False
    #server.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=True)
    #server.config.update(SESSION_COOKIE_SAMESITE="None")
    #server.config.update(SESSION_COOKIE_SECURE=True)
    atexit.register(clean_tempfiles_dir)
    CORS(server, supports_credentials=True)
    server.secret_key = 'your_secret_key'
    server.config['UPLOAD_FOLDER'] = 'pdb'
    server.config['ALLOWED_EXTENSIONS'] = {'pdb', 'pdb1'}
    server.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=1)
    server.config['TEMPFILE_DIR'] = "tempfiles"
    os.makedirs(server.config['TEMPFILE_DIR'], exist_ok=True)
    server.run(host="localhost", port=port)