import sys
import subprocess
from PyQt5.QtCore import QUrl
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtWebEngineWidgets import QWebEngineView

class NGLViewerApp(QMainWindow):
    def __init__(self):
        super(NGLViewerApp, self).__init__()

        self.setGeometry(0, 0, 1200, 900)
        self.setWindowTitle('Protein3D')

        self.web_view = QWebEngineView(self)
        self.setCentralWidget(self.web_view)

        # Enable developer tools
        dev_tools_page = QWebEngineView()
        self.web_view.page().setDevToolsPage(dev_tools_page.page())

        url = QUrl("http://127.0.0.1:8000/")
        self.web_view.load(url)
