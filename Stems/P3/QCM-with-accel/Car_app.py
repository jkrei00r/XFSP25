# region imports
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from Car_GUI import Ui_Form
from QuarterCarModel import CarController
import sys


class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        # Modify layout to add tab widget
        self.tabWidget = qtw.QTabWidget()
        self.tab_positions = qtw.QWidget()
        self.tab_forces = qtw.QWidget()

        # Clear original plot area
        while self.layout_Plot.count():
            item = self.layout_Plot.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        # Create tab structure
        self.tabWidget.addTab(self.tab_positions, "Positions")
        self.tabWidget.addTab(self.tab_forces, "Forces")
        self.layout_Plot.addWidget(self.tabWidget)

        # Set up tab layouts
        self.layout_pos = qtw.QVBoxLayout(self.tab_positions)
        self.layout_forces = qtw.QVBoxLayout(self.tab_forces)

        # Initialize controller with proper widgets
        input_widgets = (
            self.le_m1, self.le_v, self.le_k1, self.le_c1,
            self.le_m2, self.le_k2, self.le_ang, self.le_tmax,
            self.chk_IncludeAccel
        )

        display_widgets = (
            self.gv_Schematic, self.chk_LogX, self.chk_LogY,
            self.chk_LogAccel, self.chk_ShowAccel, self.lbl_MaxMinInfo,
            self.layout_pos, self.layout_forces
        )

        self.controller = CarController((input_widgets, display_widgets))

        # Connect signals
        self.btn_calculate.clicked.connect(self.controller.calculate)
        self.pb_Optimize.clicked.connect(self.doOptimize)
        self.chk_LogX.stateChanged.connect(self.controller.doPlot)
        self.chk_LogY.stateChanged.connect(self.controller.doPlot)
        self.chk_LogAccel.stateChanged.connect(self.controller.doPlot)
        self.chk_ShowAccel.stateChanged.connect(self.controller.doPlot)
        self.controller.setupEventFilter(self)
        self.controller.setupCanvasMoveEvent(self)

        self.show()

    def eventFilter(self, obj, event):
        if obj == self.gv_Schematic.scene():
            et = event.type()
            if et == qtc.QEvent.GraphicsSceneMouseMove:
                scenePos = event.scenePos()
                strScene = f"Mouse: x={scenePos.x():.1f}, y={-scenePos.y():.1f}"
                self.setWindowTitle(strScene)
            if et == qtc.QEvent.GraphicsSceneWheel:
                zm = self.controller.getZoom()
                zm += 0.1 if event.delta() > 0 else -0.1
                zm = max(0.1, zm)
                self.controller.setZoom(zm)
        self.controller.updateSchematic()
        return super().eventFilter(obj, event)

    def mouseMoveEvent_Canvas(self, event):
        if event.inaxes:
            self.controller.animate(event.xdata)
            ywheel, ybody, yroad, accel = self.controller.getPoints(event.xdata)
            self.setWindowTitle(
                f't: {event.xdata:.2f}s | '
                f'Road: {yroad * 1000:.1f}mm | '
                f'Wheel: {ywheel * 1000:.1f}mm | '
                f'Body: {ybody * 1000:.1f}mm | '
                f'Accel: {accel:.2f}g'
            )

    def doOptimize(self):
        self.controller.OptimizeSuspension()


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Quarter Car Model')
    sys.exit(app.exec())
