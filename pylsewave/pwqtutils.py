__author__ = "Georgios E. Ragkousis"
from PyQt5.QtWidgets import (QWidget, QProgressBar,
                             QPushButton, QApplication, QLabel)
from PyQt5.QtCore import QBasicTimer
from PyQt5 import QtGui
import sys
import os


class pylseWaveProgressBar(QWidget):

    def __init__(self, dt, T):
        super(pylseWaveProgressBar, self).__init__()
        self.dt = dt
        self.T = T
        self.initUI()

    def initUI(self):

        self.pbar = QProgressBar(self)
        # Creating a label
        self.progressLabel = QLabel('pylsewave solver...', self)
        self.pbar.setGeometry(30, 40, 200, 25)
        self.pbar.setMaximum(100)
        self.pbar.setMinimum(0)
        self.pbar.setRange(0, 100)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("../pylsewave.png"),
                       QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.setWindowIcon(icon)
        # self.btn = QPushButton('Start', self)
        # self.btn.move(40, 80)
        # self.btn.clicked.connect(self.doAction)

        # self.timer = QBasicTimer()
        self.t = 0
        self.step = 0
        self.setGeometry(300, 300, 280, 170)
        self.setWindowTitle('QProgressBar')
        # self.show()

    # def setValue(self, val):
    #     val = float(min(max(val, 0), 1))
    #     self._value = -270 * val
    #     self.update()

    def setProgresLabel(self):
        self.progressLabel.setText( '%0.2f%%' % self.step)

    # def timerEvent(self, e):
    #     if self.step >= self.T or self.isVisible() is not True:
    #         self.timer.stop()
    # #         self.btn.setText('Finished')
    #         return
    #
    #     self.t += self.dt
    #     self.step = (100*((self.t/self.T)))
    #     self.pbar.setFormat("%0.2f%%" % self.step)
    #     self.pbar.setValue(self.step)
        # QApplication.processEvents()

    def doAction(self):
        if self.t == 0:
            self.btn.setText('Start')
            self.t = self.dt
            while self.t <= self.T:
                self.step = (100*((self.t/self.T)))
                self.pbar.setFormat("%0.2f%%" % self.step)
                self.pbar.setValue(self.step)

                self.t += self.dt
            self.btn.setText('Finished')

        elif self.pbar.value() == 100:
            self.btn.setText('Finished')

    def closeEvent(self, e):
        # Your desired functionality here
        print('Close button pressed')
        import sys
        sys.exit(0)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = pylseWaveProgressBar(dt=1.0, T=10000.0)
    ex.show()
    ex.t = ex.dt
    while ex.t <= ex.T and ex.isVisible():
        ex.step = (100 * ((ex.t / ex.T)))
        ex.pbar.setFormat("%0.2f%%" % ex.step)
        print(ex.step)
        ex.pbar.setValue(ex.step)
        # ex.progressLabel.setText('%0.2f%%' % ex.step)
        ex.t += ex.dt
        app.processEvents()
    # ex.close()
    # QApplication.processEvents()
    print("closed")
    ret = app.exec_()
    sys.exit(ret)