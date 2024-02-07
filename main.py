from flask import Flask, render_template, request, redirect, url_for, flash, send_file
import os
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'pdb'
app.config['STRUCTURE_FILENAME'] = 'structure.pdb'
app.config['ALLOWED_EXTENSIONS'] = {'pdb', 'pdb1'}


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/pdb/structure.pdb")
def serve_pdb_file():
    return send_file("pdb/structure.pdb", as_attachment=True)


@app.route("/pdb/load_structure", methods=['POST'])
def load_structure():
    if 'file' not in request.files:
        return redirect(request.url)

    file = request.files['file']

    if file.filename == '':
        return redirect(request.url)

    if file and allowed_file(file.filename):
        filename = app.config['STRUCTURE_FILENAME']
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        return redirect(url_for('index'))
    else:
        return redirect(request.url)


if __name__ == "__main__":
    app.run(host="localhost", port=8000, debug=True)
