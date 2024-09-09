import sys
import threading
from PyQt5.QtWidgets import QApplication
from interface.app import NGLViewerApp
from server import main as run_server


def run_pyqt5_app():
    app_ = QApplication(sys.argv)
    viewer_app = NGLViewerApp()
    viewer_app.show()
    sys.exit(app_.exec_())


def run_flask_server():
    run_server.main()

pyqt5_thread = threading.Thread(target=run_pyqt5_app)
flask_thread = threading.Thread(target=run_flask_server)

pyqt5_thread.start()
flask_thread.setDaemon(True)
flask_thread.start()