import sys
import subprocess
from PyQt5.QtCore import QUrl
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEnginePage, QWebEngineDownloadItem

class NGLViewerApp(QMainWindow):
    def __init__(self):
        super(NGLViewerApp, self).__init__()

        self.setWindowTitle('Protein3D')
        self.web_view = QWebEngineView(self)
        self.page = QWebEnginePage(self)
        self.web_view.setPage(self.page)
        self.setCentralWidget(self.web_view)
        url = QUrl("http://127.0.0.1:8000/")
        self.web_view.load(url)
        self.showMaximized()
