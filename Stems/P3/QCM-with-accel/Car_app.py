# region imports
import sys
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from Car_GUI import Ui_Form
from QuarterCarModel import CarController


class MainWindow(qtw.QWidget, Ui_Form):
    """Main application window handling GUI setup and controller integration"""

    def __init__(self):
        super().__init__()
        try:
            # Initialize UI components from generated code
            self.setupUi(self)

            # Configure tab system for plots
            self._setup_tabs()

            # Initialize controller with UI components
            self.controller = self._create_controller()

            # Connect UI signals to controller methods
            self._connect_signals()

            self.show()  # Display the main window
            print("GUI initialized successfully")

        except Exception as e:
            qtw.QMessageBox.critical(None, "Startup Error", f"Failed to initialize:\n{str(e)}")
            sys.exit(1)

    def _setup_tabs(self):
        """Configure tab widget for position/force plots"""
        # Create tab container and individual tabs
        self.tabWidget = qtw.QTabWidget()
        self.tab_positions = qtw.QWidget()
        self.tab_forces = qtw.QWidget()

        # Clear existing items in plot area
        while self.layout_Plot.count():
            item = self.layout_Plot.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        # Configure tab structure
        self.tabWidget.addTab(self.tab_positions, "Positions")
        self.tabWidget.addTab(self.tab_forces, "Forces")
        self.layout_Plot.addWidget(self.tabWidget)

        # Create layouts for each tab
        self.layout_pos = qtw.QVBoxLayout(self.tab_positions)
        self.layout_forces = qtw.QVBoxLayout(self.tab_forces)

    def _create_controller(self):
        """Initialize controller with UI components"""
        input_widgets = (
            self.le_m1,  # 0: Car body mass
            self.le_v,  # 1: Vehicle speed
            self.le_k1,  # 2: Suspension spring
            self.le_c1,  # 3: Damper coefficient
            self.le_m2,  # 4: Wheel mass
            self.le_k2,  # 5: Tire spring
            self.le_ang,  # 6: Ramp angle
            self.le_tmax,  # 7: Max simulation time
            self.chk_IncludeAccel  # 8: Acceleration checkbox
        )

        display_widgets = (
            self.gv_Schematic,  # Schematic diagram view
            self.chk_LogX,  # X-axis log scale
            self.chk_LogY,  # Y-axis log scale
            self.chk_LogAccel,  # Acceleration log scale
            self.chk_ShowAccel,  # Show acceleration
            self.lbl_MaxMinInfo,  # Min/max display label
            self.layout_pos,  # Position plot layout
            self.layout_forces  # Force plot layout
        )

        return CarController((input_widgets, display_widgets))

    def _connect_signals(self):
        """Connect UI element signals to controller methods"""
        self.btn_calculate.clicked.connect(self.controller.calculate)
        self.pb_Optimize.clicked.connect(self.controller.OptimizeSuspension)
        # Connect plot configuration checkboxes
        self.chk_LogX.stateChanged.connect(self.controller.doPlot)
        self.chk_LogY.stateChanged.connect(self.controller.doPlot)
        self.chk_LogAccel.stateChanged.connect(self.controller.doPlot)
        self.chk_ShowAccel.stateChanged.connect(self.controller.doPlot)

    def eventFilter(self, obj, event):
        """Handle mouse events for schematic view interaction"""
        # ... (existing event filter code) ...


if __name__ == '__main__':
    # Configure application entry point
    app = qtw.QApplication(sys.argv)
    try:
        window = MainWindow()
        window.setWindowTitle('Quarter Car Model')
        sys.exit(app.exec())
    except Exception as e:
        qtw.QMessageBox.critical(None, "Fatal Error",
                                 f"Application failed to start:\n{str(e)}")
        sys.exit(1)
