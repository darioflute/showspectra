#!/usr/bin/env python

import sys,os
import numpy as np


class GUI (QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'Show Spectra'
        self.left = 10
        self.top = 10
        self.width = 640
        self.height = 480

        # Ordinate closing
        self.setAttribute(Qt.WA_DeleteOnClose)

        # Path of the package
        self.path0, file0 = os.path.split(__file__)

        # Style

        self.initUI()

    def initUI(self):
        """ User interface """

        self.setWindowTitle(self.title)
        self.setGeometry(self.left,self.top,self.width,self.height)

        # Main widget
        wid = QWidget()
        self.setCentralWidget(wid)

        # Main layout
        mainLayout = QHBoxLayout()
        
        # main panel
        self.createImagePanel()

        # Add panel to main layout
        mainLayout.addWidget(self.imagePanel)
        wid.setLayout(mainLayout)
        self.show()

        # Menu
        self.createMenu()

        
