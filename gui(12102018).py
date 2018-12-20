# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gui.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!
import os 

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import (QApplication, QDialog,
                             QProgressBar, QPushButton, QMessageBox)
import matplotlib.pyplot as plt
from matplotlib import style
import T2H, PLOT
import flopy
from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
else:
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure

#%%
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("T2H Graphical User Interface")
        MainWindow.resize(1280, 800)
        self.centralWidget = QtWidgets.QWidget(MainWindow)
        self.centralWidget.setObjectName("centralWidget")

#%% QFrames
        self.frame_1 = QtWidgets.QFrame(self.centralWidget)
        self.frame_1.setGeometry(QtCore.QRect(810, 70, 461, 201))
        self.frame_1.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_1.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_1.setObjectName("frame_2")
        
        self.frame_2 = QtWidgets.QFrame(self.centralWidget)
        self.frame_2.setGeometry(QtCore.QRect(810, 280, 461, 101))
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")

        self.frame_3 = QtWidgets.QFrame(self.centralWidget)
        self.frame_3.setGeometry(QtCore.QRect(810, 390, 461, 31))
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")                
#%% QLabels
        self.sedK = QtWidgets.QLabel(self.frame_2)
        self.sedK.setGeometry(QtCore.QRect(30, 10, 141, 16))
        self.sedK.setObjectName("sedK")
        
        self.aqK = QtWidgets.QLabel(self.frame_2)
        self.aqK.setGeometry(QtCore.QRect(30, 40, 141, 16))
        self.aqK.setObjectName("aqK")
        
        self.faultK = QtWidgets.QLabel(self.frame_2)
        self.faultK.setGeometry(QtCore.QRect(30, 70, 141, 16))
        self.faultK.setObjectName("faultK")   

        self.sedKN = QtWidgets.QLabel(self.centralWidget)
        self.sedKN.setGeometry(QtCore.QRect(910, 500, 141, 16))
        self.sedKN.setObjectName("sedKN")
        
        self.sedKNlabel = QtWidgets.QLabel(self.centralWidget)
        self.sedKNlabel.setGeometry(QtCore.QRect(1100, 500, 61, 16))
        self.sedKNlabel.setObjectName("sedKNlabel")
        
        self.aquiferKNlabel = QtWidgets.QLabel(self.centralWidget)
        self.aquiferKNlabel.setGeometry(QtCore.QRect(1100, 520, 61, 16))
        self.aquiferKNlabel.setObjectName("aquiferKNlabel")
        
        self.aqKN = QtWidgets.QLabel(self.centralWidget)
        self.aqKN.setGeometry(QtCore.QRect(910, 520, 81, 16))
        self.aqKN.setObjectName("aqKN")
        
        self.faultKN = QtWidgets.QLabel(self.centralWidget)
        self.faultKN.setGeometry(QtCore.QRect(910, 540, 81, 16))
        self.faultKN.setObjectName("faultKN")
        
        self.faultKNlabel = QtWidgets.QLabel(self.centralWidget)
        self.faultKNlabel.setGeometry(QtCore.QRect(1100, 540, 61, 16))
        self.faultKNlabel.setObjectName("faultKNlabel")

        self.label_21 = QtWidgets.QLabel(self.frame_3)
        self.label_21.setGeometry(QtCore.QRect(10, 7, 141, 16))
        self.label_21.setObjectName("label_21")

        self.visoptionsLabel = QtWidgets.QLabel(self.centralWidget)
        self.visoptionsLabel.setGeometry(QtCore.QRect(20, 540, 141, 16))
        self.visoptionsLabel.setObjectName("visoptionsLabel")
        
        self.fileLabel = QtWidgets.QLabel(self.centralWidget)
        self.fileLabel.setGeometry(QtCore.QRect(810, 4, 60, 16))
        self.fileLabel.setObjectName("fileLabel")
        
        self.fileLabel_path = QtWidgets.QLabel(self.centralWidget)
        self.fileLabel_path.setGeometry(QtCore.QRect(880, 4, 320, 16))
        self.fileLabel_path.setObjectName("fileLabel_path")
        
        self.label = QtWidgets.QLabel(self.centralWidget)
        self.label.setGeometry(QtCore.QRect(814, 51, 241, 16))
        self.label.setObjectName("label")
        
        self.nz = QtWidgets.QLabel(self.centralWidget)
        self.nz.setGeometry(QtCore.QRect(840, 104, 141, 16))
        self.nz.setObjectName("nz")

        self.targetperiod = QtWidgets.QLabel(self.centralWidget)
        self.targetperiod.setGeometry(QtCore.QRect(840, 80, 151, 16))
        self.targetperiod.setObjectName("targetperiod")
        
        self.nzfixed = QtWidgets.QLabel(self.centralWidget)
        self.nzfixed.setGeometry(QtCore.QRect(840, 128, 141, 16))
        self.nzfixed.setObjectName("nzfixed")

        self.constrecharge = QtWidgets.QLabel(self.centralWidget)
        self.constrecharge.setGeometry(QtCore.QRect(840, 176, 151, 16))
        self.constrecharge.setObjectName("constrecharge")

        #
        self.hiniratio = QtWidgets.QLabel(self.centralWidget)
        self.hiniratio.setGeometry(QtCore.QRect(840, 242, 151, 16))
        self.hiniratio.setObjectName("hiniratio")

        self.datvar = QtWidgets.QLabel(self.centralWidget)
        self.datvar.setGeometry(QtCore.QRect(840, 152, 161, 16))
        self.datvar.setObjectName("datvar")
        

        # Recharge input 
        self.constrecharge_2 = QtWidgets.QLabel(self.centralWidget)
        self.constrecharge_2.setGeometry(QtCore.QRect(840, 200, 151, 16))
        self.constrecharge_2.setObjectName("constrecharge_2")
        
        # Image pane
        self.image = QtWidgets.QLabel(self.centralWidget)
        self.image.setGeometry(QtCore.QRect(10, 10, 780, 520))
        self.image.setObjectName("image")
        self.pixmap = QtGui.QPixmap("logo.png")
        self.image.setPixmap(self.pixmap)
#%% QLineEdits
        self.sedKlineEdit = QtWidgets.QLineEdit(self.frame_2)
        self.sedKlineEdit.setGeometry(QtCore.QRect(260, 10, 113, 21))
        self.sedKlineEdit.setObjectName("sedKlineEdit")
        self.sedKlineEdit.setText("547.5")
        #
        self.aqKlineEdit = QtWidgets.QLineEdit(self.frame_2)
        self.aqKlineEdit.setGeometry(QtCore.QRect(260, 40, 113, 21))
        self.aqKlineEdit.setObjectName("aqKlineEdit")
        self.aqKlineEdit.setText("36.5")
        #
        self.faultKlineEdit = QtWidgets.QLineEdit(self.frame_2)
        self.faultKlineEdit.setGeometry(QtCore.QRect(260, 70, 113, 21))
        self.faultKlineEdit.setObjectName("faultKlineEdit")
        self.faultKlineEdit.setText("0.0365")
        #
        self.nzfline = QtWidgets.QLineEdit(self.centralWidget)
        self.nzfline.setGeometry(QtCore.QRect(1070, 128, 113, 21))
        self.nzfline.setObjectName("nzfline")
        self.nzfline.setText("10")
        #
        self.nzline = QtWidgets.QLineEdit(self.centralWidget)
        self.nzline.setGeometry(QtCore.QRect(1070, 104, 113, 21))
        self.nzline.setObjectName("nzline")
        self.nzline.setText("40")
        #
        self.datline = QtWidgets.QLineEdit(self.centralWidget)
        self.datline.setGeometry(QtCore.QRect(1070, 152, 113, 21))
        self.datline.setObjectName("datline")
        self.datline.setText("-10000")
        #
        self.hiniratioLineEdit = QtWidgets.QLineEdit(self.centralWidget)
        self.hiniratioLineEdit.setGeometry(QtCore.QRect(1070, 242, 113, 21))
        self.hiniratioLineEdit.setObjectName("hiniratioLineEdit")
        self.hiniratioLineEdit.setText("0.9")
        
        #
        self.datvarline = QtWidgets.QLineEdit(self.centralWidget)
        self.datvarline.setGeometry(QtCore.QRect(1070, 176, 113, 21))
        self.datvarline.setObjectName("datvarline")
        self.datvarline.setText("-3000")

        self.rchline = QtWidgets.QLineEdit(self.centralWidget)
        self.rchline.setGeometry(QtCore.QRect(1070, 200, 113, 21))
        self.rchline.setObjectName("rchline")
        self.rchline.setText("0.05")
        
        # Ma input lineedit
        self.maline = QtWidgets.QLineEdit(self.centralWidget)
        self.maline.setGeometry(QtCore.QRect(1070, 80, 113, 21))
        self.maline.setObjectName("maline")
        self.maline.setText("12.5")
        
#%% QPushButtons
        self.load = QtWidgets.QPushButton(self.centralWidget)
        self.load.setGeometry(QtCore.QRect(1100, -1, 71, 32))
        self.load.setObjectName("loadButton")
        self.load.clicked.connect(self.fileloader)   
        
        self.load1 = QtWidgets.QPushButton(self.centralWidget)
        self.load1.setGeometry(QtCore.QRect(1170, -1, 101, 32))
        self.load1.setObjectName("loadButton1")
        self.load1.clicked.connect(self.fileloader)   
        
        self.applyButton = QtWidgets.QPushButton(self.frame_1)
        self.applyButton.setGeometry(QtCore.QRect(380, 60, 81, 81))
        self.applyButton.setObjectName("applyButton")
        self.applyButton.clicked.connect(self.applyclicked)   

        self.fileDialog_3 = QtWidgets.QPushButton(self.frame_2)
        self.fileDialog_3.setGeometry(QtCore.QRect(380, 20, 81, 71))
        self.fileDialog_3.setObjectName("fileDialog_3")
        self.fileDialog_3.clicked.connect(self.applyCalClicked) 

        # Model run button
        self.ModelRunButton = QtWidgets.QPushButton(self.centralWidget)
        self.ModelRunButton.setGeometry(QtCore.QRect(640, 620, 113, 32))
        self.ModelRunButton.setObjectName("ModelRunButton")
        self.ModelRunButton.clicked.connect(self.run)
        
        self.QuitButton = QtWidgets.QPushButton(self.centralWidget)
        self.QuitButton.setGeometry(QtCore.QRect(760, 620, 113, 32))
        self.QuitButton.setObjectName("QuitButton")
        self.QuitButton.clicked.connect(QCoreApplication.instance().quit)
        
        self.VtkOutputButton = QtWidgets.QPushButton(self.centralWidget)
        self.VtkOutputButton.setGeometry(QtCore.QRect(880, 620, 113, 32))
        self.VtkOutputButton.setObjectName("VtkOutputButton")
#        self.VtkOutputButton.clicked.connect(self.vtk)
        
        self.PlotButton = QtWidgets.QPushButton(self.centralWidget)
        self.PlotButton.setGeometry(QtCore.QRect(460, 560, 113, 32))
        self.PlotButton.setObjectName("PlotButton")
        self.PlotButton.clicked.connect(self.plot)
#%% QGraphicsViews    
        self.figure = plt.figure(figsize=(12,12))
        self.canvas = FigureCanvas(self.figure)

#%% QComboBoxes
        # File combo box
        self.fileBox = QtWidgets.QComboBox(self.centralWidget)
        self.fileBox.setGeometry(QtCore.QRect(808, 25, 461, 26))
        self.fileBox.setObjectName("fileBox")
        # Solver selection combo box
        self.solverBox = QtWidgets.QComboBox(self.frame_3)
        self.solverBox.setGeometry(QtCore.QRect(63, 2, 281, 26))
        self.solverBox.setObjectName("solverBox")
        self.solverBox.addItem("xMD")
        self.solverBox.addItem("GMRES")
        #
        self.visComboBox = QtWidgets.QComboBox(self.centralWidget)
        self.visComboBox.setGeometry(QtCore.QRect(10, 560, 441, 26))
        self.visComboBox.setObjectName("visComboBox")
        self.visComboBox.addItem("Cross Section")
        self.visComboBox.addItem("Fault Plane")
        self.visComboBox.addItem("Vertical Flow Barriers (VFB)")
        self.visComboBox.addItem("Horizontal Flow Barriers (HFB)")
        
#%% QCheckBoxes
        #
        self.elevdependentChecker = QtWidgets.QCheckBox(self.centralWidget)
        self.elevdependentChecker.setGeometry(QtCore.QRect(860, 220, 231, 20))
        self.elevdependentChecker.setObjectName("elevdependentChecker")

#%% QProgressBars
        self.progress = QProgressBar(self.centralWidget)
        self.progress.setGeometry(10, 620, 600, 25)   
        self.progress.setMaximum(100)

#%% Mainwindows       
        MainWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QtWidgets.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 1024, 22))
        self.menuBar.setObjectName("menuBar")
        self.menuT2H_Main = QtWidgets.QMenu(self.menuBar)
        self.menuT2H_Main.setObjectName("menuT2H_Main")
        self.menuT2H_Checker = QtWidgets.QMenu(self.menuBar)
        self.menuT2H_Checker.setObjectName("menuT2H_Checker")
        self.menuT2H_Plot = QtWidgets.QMenu(self.menuBar)
        self.menuT2H_Plot.setObjectName("menuT2H_Plot")
        
        MainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QtWidgets.QToolBar(MainWindow)
        self.mainToolBar.setObjectName("mainToolBar")
        
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtWidgets.QStatusBar(MainWindow)
        self.statusBar.setObjectName("statusBar")
        
        MainWindow.setStatusBar(self.statusBar)
        self.menuBar.addAction(self.menuT2H_Main.menuAction())
        self.menuBar.addAction(self.menuT2H_Checker.menuAction())
        self.menuBar.addAction(self.menuT2H_Plot.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
#%% Functions    
    def applyclicked(self):
        self.Ma = float(self.maline.text())
        self.Ma = format(self.Ma, '.1f')
        self.nz = int(self.nzline.text())
        self.nz_fixed = int(self.nzfline.text())
        self.dx = 1000
        self.dy = 1000
        self.inz = self.nz - self.nz_fixed
        self.dat = int(self.datline.text())
        self.dat_var = int(self.datvarline.text())
        self.idat = self.dat - self.dat_var
        self.rech = float(self.rchline.text())
        self.perm_sed = float(self.sedKlineEdit.text())
        self.hratio = float(self.hiniratioLineEdit.text())
        self.Kconst = float(self.aqKlineEdit.text())
        self.hydchr = self.Kconst/1000
        self.target_row = 101
        self.iskip = 4
        self.ivtk = 1
        self.h_tol = 1e-4
        self.fileLabel_path.setText("/tisc_output/topo_" + self.Ma +"0Ma.txt")
        self.ans = QMessageBox.question(self.centralWidget, "Confirmation",\
                          "Are these correct?\n" + "Period: " + self.Ma\
                          + "Ma\n" + "Nz: " + str(self.nz) +"\n" + "Datum: "\
                          + str(self.dat) + " m\n", QMessageBox.Yes, QMessageBox.No)
        if self.ans == QMessageBox.Yes:
            self.rchline.setEnabled(False)
            self.maline.setEnabled(False)
            self.nzline.setEnabled(False)
            self.nzfline.setEnabled(False)
            self.datline.setEnabled(False)
            self.datvarline.setEnabled(False)
            self.hiniratioLineEdit.setEnabled(False)
            QMessageBox.about(self.centralWidget, "Confirmed", "Properties confirmed")
        else:
            QMessageBox.about(self.centralWidget, "Check values", "Check values again!")

    def applyCalClicked(self):
        self.perm_sed = self.sedKlineEdit.text()
        self.Kconst = self.aqKlineEdit.text()
        self.hydchr = self.faultKlineEdit.text()
        self.sedKNlabel.setText(str(float(self.perm_sed)/float(self.rchline.text())))
        self.aquiferKNlabel.setText(str(float(self.Kconst)/float(self.rchline.text())))
        self.faultKNlabel.setText(str(float(self.hydchr)/float(self.rchline.text())))
        self.ans = QMessageBox.question(self.centralWidget, "Confirmation",\
                          "Are these correct?\n" + "Period: " + self.Ma\
                          + "Ma\n" + "Nz: " + str(self.nz) +"\n" + "Datum: "\
                          + str(self.dat) + " m\n", QMessageBox.Yes, QMessageBox.No)
        if self.ans == QMessageBox.Yes:
            self.sedKlineEdit.setEnabled(False)
            self.aqKlineEdit.setEnabled(False)
            self.faultKlineEdit.setEnabled(False)
            QMessageBox.about(self.centralWidget, "Confirmed", "Properties confirmed")
        else:
            QMessageBox.about(self.centralWidget, "Check values", "Check values again!")
        
#%%    
    def run(self):
        self.Ma = float(self.maline.text())
        self.Ma = format(self.Ma, '.1f')
        self.nz = int(self.nzline.text())
        self.nz_fixed = int(self.nzfline.text())
        self.dx = 1000
        self.dy = 1000
        self.inz = self.nz - self.nz_fixed
        self.dat = int(self.datline.text())
        self.dat_var = int(self.datvarline.text())
        self.idat = self.dat - self.dat_var
        self.rech = float(self.rchline.text())
        self.perm_sed = float(self.sedKlineEdit.text())
        self.hratio = float(self.hiniratioLineEdit.text())
        self.Kconst = float(self.aqKlineEdit.text())
        self.hydchr = self.Kconst/1000
        self.target_row = 101
        self.iskip = 4
        self.ivtk = 1
        self.h_tol = 1e-4
        self.model = T2H.main(self.Ma, self.nz, self.nz_fixed, self.inz, self.dx,\
                           self.dy, self.dat, self.dat_var, self.idat\
                           , self.rech, self.perm_sed, self.target_row,\
                           self.Kconst, self.hratio, self.hydchr,\
                           self.iskip, self.ivtk, self.h_tol)
        self.mf = self.model.mf
        self.mf.dis.check()
        self.mf.write_input()
        self.mf.run_model()
        return self.mf
        
    def plot(self):
        try: 
            self.mf
        except AttributeError:
            QMessageBox.about(self.centralWidget, "Warning", "Please run a model first")
        else:
            self.vcb = self.visComboBox.itemData
            print(self.vcb)
            if self.vcb == "Cross Section":
                figheadxsect, axheadxsect = plt.subplots(figsize=(40,5))
                self.mfxsect = PLOT.fmfxsect(self.mf, self.model.mfdis, self.target_row, axheadxsect).mfxsect
                self.a = PLOT.head(self.mf, self.model.fdirmodel).a
                self.headc = PLOT.headc(self.mfxsect, self.a)
                self.headcontour = self.headc.headcontour
                self.gdplot = self.mfxsect.plot_grid(color='r', linewidths=0.2)
                self.BCplot = self.mfxsect.plot_ibound(self.model.ibound, color_noflow = 'black',\
                                                       color_ch = 'blue', head = self.a)
                self.canvas.draw()
            
                print("plot")
    
    def fileloader(self):
        self.path = os.getcwd() + "/tisc_output/"
        self.l = os.listdir(self.path)
        self.bdtopo = [0]*len(self.l)
        self.topo = [0]*len(self.l)
        self.fault = [0]*len(self.l)
        self.sedthick = [0]*len(self.l)
        for file in range(len(self.l)):
            if self.l[file].startswith("bdtopo"):
                if os.stat(self.path+self.l[file]).st_size > 5: # greater than 5 bytes
                    self.bdtopo[file] = float(self.l[file][7:]\
                               .split("Ma.txt")[0])
            elif self.l[file].startswith("topo"):
                if os.stat(self.path+self.l[file]).st_size > 5: # greater than 5 bytes
                    self.topo[file] = float(self.l[file][5:]\
                             .split("Ma.txt")[0])
            elif self.l[file].startswith("fault"):
                if os.stat(self.path+self.l[file]).st_size > 5: # greater than 5 bytes
                    self.fault[file] = float(self.l[file][6:]\
                              .split("Ma.txt")[0])
            elif self.l[file].startswith("sedthick"):
                if os.stat(self.path+self.l[file]).st_size > 5: # greater than 5 bytes
                    self.sedthick[file] = float(self.l[file][9:]\
                            .split("Ma.txt")[0])
        self.a = list(filter((0).__ne__, self.topo))
        self.a.sort()
        self.b = list(filter((0).__ne__, self.bdtopo))
        self.b.sort()
        self.c = list(filter((0).__ne__, self.fault))
        self.c.sort()
        self.d = list(filter((0).__ne__, self.sedthick))
        self.d.sort()
        self.df = []
        for nfile in range(len(self.a)):
            if self.b.count(self.a[nfile]) == 1:
                if self.c.count(self.a[nfile]) == 1:
                    if self.d.count(self.a[nfile]) == 1:
                        data = [self.a[nfile], "y", "y", "y", "y"]
                        self.df.append(data)
                    elif self.d.count(self.a[nfile]) == 0:
                        data = [self.a[nfile], "y", "y", "y", "n"]
                        self.df.append(data)
                elif self.c.count(self.a[nfile]) == 0:
                    if self.d.count(self.a[nfile]) == 1:
                        data = [self.a[nfile], "y", "y", "n", "y"]
                        self.df.append(data)
                    elif self.d.count(self.a[nfile]) == 0:
                        data = [self.a[nfile], "y", "y", "n", "n"]
                        self.df.append(data)                    
            elif self.b.count(self.a[nfile]) == 0:
                if self.c.count(self.a[nfile]) == 1:
                    if self.d.count(self.a[nfile]) == 1:
                        data = [self.a[nfile], "y", "n", "y", "y"]
                        self.df.append(data)
                    elif self.d.count(self.a[nfile]) == 0:
                        data = [self.a[nfile], "y", "n", "y", "n"]
                        self.df.append(data)
                elif self.c.count(self.a[nfile]) == 0:
                    if self.d.count(self.a[nfile]) == 1:
                        data = [self.a[nfile], "y", "n", "n", "y"]
                        self.df.append(data)
                    elif self.d.count(self.a[nfile]) == 0:
                        data = [self.a[nfile], "y", "n", "n", "n"]
                        self.df.append(data)
        for age in range(len(self.a)):
            if self.df[age][2] == "y" and self.df[age][3] == "y" and self.df[age][4] == "y":
                self.fileBox.addItem("Snapshot:" + str(self.df[age][0]) + "Ma | Faults | Sediments")
            elif self.df[age][2] == "y" and self.df[age][3] == "y" and self.df[age][4] == "n":
                self.fileBox.addItem("Snapshot:" + str(self.df[age][0]) + "Ma | Faults | No Sediments")
            elif self.df[age][2] == "y" and self.df[age][3] == "n" and self.df[age][4] == "y":
                self.fileBox.addItem("Snapshot:" + str(self.df[age][0]) + "Ma | No Faults | Sediments")
            elif self.df[age][2] == "y" and self.df[age][3] == "n" and self.df[age][4] == "n":
                self.fileBox.addItem("Snapshot:" + str(self.df[age][0]) + "Ma | No Faults | No Sediments")

 
#%%        

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "T2H Graphical User Interface"))
        self.applyButton.setText(_translate("MainWindow", "Apply"))
        self.sedK.setText(_translate("MainWindow", "Sediment K (m/yr)"))
        self.aqK.setText(_translate("MainWindow", "Aquifer K (m/yr)"))
        self.faultK.setText(_translate("MainWindow", "Fault zone K (m/yr)"))
        self.fileDialog_3.setText(_translate("MainWindow", "Apply"))
        self.sedKN.setText(_translate("MainWindow", "Sediment K / N:"))
        self.sedKNlabel.setText(_translate("MainWindow", "N/A"))
        self.aquiferKNlabel.setText(_translate("MainWindow", "N/A"))
        self.aqKN.setText(_translate("MainWindow", "Aquifer K / N:"))
        self.faultKN.setText(_translate("MainWindow", "Fault K / N:"))
        self.faultKNlabel.setText(_translate("MainWindow", "N/A"))
        self.label_21.setText(_translate("MainWindow", "Solver"))
        self.ModelRunButton.setText(_translate("MainWindow", "Execute"))
        self.load.setText(_translate("MainWindow", "Load"))
        self.load1.setText(_translate("MainWindow", "Set selected"))
        self.QuitButton.setText(_translate("MainWindow", "Abort"))
        self.VtkOutputButton.setText(_translate("MainWindow", "VTK output"))
        self.PlotButton.setText(_translate("MainWindow", "Plot"))
        self.visoptionsLabel.setText(_translate("MainWindow", "Visualization options"))
        self.fileLabel.setText(_translate("MainWindow", "File: "))
        self.fileLabel_path.setText(_translate("MainWindow", "path"))
        self.label.setText(_translate("MainWindow", "*dx = dy = 1,000 m fixed in this version"))
        self.nz.setText(_translate("MainWindow", "Number of layers (nz)"))
        self.targetperiod.setText(_translate("MainWindow", "Target period (Ma)"))
        self.nzfixed.setText(_translate("MainWindow", "Fixed layers (nz_fixed)"))
        self.constrecharge.setText(_translate("MainWindow", "Datum of variable dz (m)"))
        self.hiniratio.setText(_translate("MainWindow", "Initial head ratio to topo."))
        self.elevdependentChecker.setText(_translate("MainWindow", "Elevation-dependent recharge"))
        self.datvar.setText(_translate("MainWindow", "Model datum (m)"))
        self.constrecharge_2.setText(_translate("MainWindow", "Const. Recharge (m/yr)"))
        self.menuT2H_Main.setTitle(_translate("MainWindow", "T2H Main"))
        self.menuT2H_Checker.setTitle(_translate("MainWindow", "T2H Checker"))
        self.menuT2H_Plot.setTitle(_translate("MainWindow", "T2H Plot"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

