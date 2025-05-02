# region imports
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from PyQt5 import QtWidgets as qtw, QtCore as qtc, QtGui as qtg
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure


class CarModel:
    """Data model storing simulation parameters and results"""

    def __init__(self):
        # Simulation parameters with default values
        self.m1 = 450  # Car body mass (kg)
        self.m2 = 20  # Wheel mass (kg)
        self.k1 = 15000  # Suspension spring constant (N/m)
        self.k2 = 90000  # Tire spring constant (N/m)
        self.c1 = 4500  # Damper coefficient (N·s/m)
        self.v = 120  # Vehicle speed (kph)
        self.yangdeg = 45  # Ramp angle (degrees)
        self.tmax = 3  # Maximum simulation time (s)
        self.accelLim = 1.5  # Acceleration limit (g)

        # Calculated parameters
        self.ymag = 0.1515  # Ramp height (m)
        self.tramp = 0  # Time to traverse ramp (s)
        self.angrad = 0  # Ramp angle in radians
        self.SSE = 0.0  # Sum of squared errors
        self.accelMax = 0  # Maximum acceleration

        # Solution arrays
        self.timeData = np.array([])
        self.results = None
        self.springForceData = []
        self.dashpotForceData = []
        self.tireForceData = []
        self.bodyPosData = []
        self.wheelPosData = []
        self.roadPosData = []
        self.accelBodyData = []


class CarView:
    """Handles visual representation of model data"""

    def __init__(self, args):
        self.input_widgets, self.display_widgets = args
        self._parse_display_widgets()
        self._initialize_plots()
        self.buildScene()
        self.zoom = 1.0

    def _parse_display_widgets(self):
        """Extract display components from passed arguments"""
        (self.gv_Schematic, self.chk_LogX, self.chk_LogY,
         self.chk_LogAccel, self.chk_ShowAccel, self.lbl_MaxMinInfo,
         self.layout_pos, self.layout_forces) = self.display_widgets

    def _initialize_plots(self):
        """Set up matplotlib figures and canvases for plotting"""
        # Position plot configuration
        self.figure_pos = Figure()
        self.canvas_pos = FigureCanvasQTAgg(self.figure_pos)
        self.ax_pos = self.figure_pos.add_subplot(111)
        self.layout_pos.addWidget(NavigationToolbar2QT(self.canvas_pos, self.canvas_pos))
        self.layout_pos.addWidget(self.canvas_pos)

        # Force plot configuration
        self.figure_force = Figure()
        self.canvas_force = FigureCanvasQTAgg(self.figure_force)
        self.ax_force = self.figure_force.add_subplot(111)
        self.layout_forces.addWidget(NavigationToolbar2QT(self.canvas_force, self.canvas_force))
        self.layout_forces.addWidget(self.canvas_force)

    def buildScene(self):
        """Initialize schematic visualization components"""
        self.scene = qtw.QGraphicsScene()
        self.scene.setSceneRect(-200, -200, 400, 400)
        self.gv_Schematic.setScene(self.scene)

        # Draw schematic components
        pen = qtg.QPen(qtc.Qt.black)
        brush = qtg.QBrush(qtc.Qt.SolidPattern)

        # Road surface
        self.scene.addLine(-150, 50, 150, 50, pen)
        # Wheel
        self.scene.addEllipse(-25, -25, 50, 50, pen, brush)
        # Suspension
        self.scene.addLine(0, -25, 0, -100, pen)
        # Car body
        self.scene.addRect(-50, -150, 100, 50, pen, brush)

    def updateView(self, model=None):
        """Update all visual elements with current model data"""
        self.doPlot(model)
        self.doForcePlot(model)
        self.updateSchematic()

    def doPlot(self, model=None):
        """Update position vs time plot"""
        self.ax_pos.clear()
        if model and model.results:
            t = model.timeData
            self.ax_pos.plot(t, model.bodyPosData, 'b-', label='Body Position')
            self.ax_pos.plot(t, model.wheelPosData, 'r-', label='Wheel Position')
            self.ax_pos.plot(t, model.roadPosData, 'k--', label='Road Profile')

            if model.accelBodyData.any() and self.chk_ShowAccel.isChecked():
                ax_accel = self.ax_pos.twinx()
                ax_accel.plot(t, model.accelBodyData, 'g-', label='Acceleration')
                ax_accel.set_ylabel('Acceleration (g)', color='g')

            self.ax_pos.set_xlabel('Time (s)')
            self.ax_pos.set_ylabel('Position (m)')
            self.ax_pos.legend()
            self.ax_pos.grid(True)
            self.canvas_pos.draw()

    def doForcePlot(self, model=None):
        """Update force vs time plot"""
        self.ax_force.clear()
        if model and model.springForceData:
            t = model.timeData
            self.ax_force.plot(t, model.springForceData, 'b-', label='Spring Force')
            self.ax_force.plot(t, model.dashpotForceData, 'r-', label='Damper Force')
            self.ax_force.plot(t, model.tireForceData, 'g-', label='Tire Force')
            self.ax_force.set_xlabel('Time (s)')
            self.ax_force.set_ylabel('Force (N)')
            self.ax_force.legend()
            self.ax_force.grid(True)
            self.canvas_force.draw()

    def updateSchematic(self):
        """Update schematic visualization based on current model state"""
        # Implementation for animating schematic view
        pass

    def animate(self, t):
        """Update schematic visualization at specific time t"""
        # Implementation for time-based animation
        pass

    def getZoom(self):
        return self.zoom

    def setZoom(self, val):
        self.zoom = max(0.1, min(2.0, val))
        self.updateSchematic()




class CarController:
    """Mediates between model and view, handles calculations and optimization"""

    def __init__(self, args):
        self.model = CarModel()
        self.view = CarView(args)
        self.chk_IncludeAccel = None

    def ode_system(self, t, y):
        """System of ODEs defining the quarter car model dynamics"""
        x1, x1_dot, x2, x2_dot = y

        # Calculate road profile
        if t < self.model.tramp:
            y_road = t * self.model.ymag / self.model.tramp
        else:
            y_road = self.model.ymag

        # Calculate forces
        F_spring = self.model.k1 * (x2 - x1)
        F_damper = self.model.c1 * (x2_dot - x1_dot)
        F_tire = self.model.k2 * (y_road - x2)

        return [
            x1_dot,
            (F_spring + F_damper) / self.model.m1,
            x2_dot,
            (-F_spring - F_damper + F_tire) / self.model.m2
        ]

    def calculate(self):
        """Main calculation routine triggered by Calculate button"""
        try:
            self._load_parameters()
            self._solve_ode()
            self._postprocess_results()
            self.view.updateView(self.model)
        except Exception as e:
            qtw.QMessageBox.critical(None, "Calculation Error", str(e))

    def _load_parameters(self):
        """Load and validate parameters from GUI inputs"""
        self.model.m1 = float(self.view.input_widgets[0].text())
        self.model.v = float(self.view.input_widgets[1].text())
        self.model.k1 = float(self.view.input_widgets[2].text())
        self.model.c1 = float(self.view.input_widgets[3].text())
        self.model.m2 = float(self.view.input_widgets[4].text())
        self.model.k2 = float(self.view.input_widgets[5].text())
        self.model.yangdeg = float(self.view.input_widgets[6].text())
        self.model.tmax = float(self.view.input_widgets[7].text())

        # Derived parameters
        v_mps = self.model.v * 1000 / 3600
        self.model.angrad = math.radians(self.model.yangdeg)
        self.model.tramp = self.model.ymag / (math.sin(self.model.angrad) * v_mps)
        self.model.timeData = np.linspace(0, self.model.tmax, 2000)

    def _solve_ode(self):
        """Solve the system of ODEs using LSODA method"""
        ic = [0.0, 0.0, 0.0, 0.0]
        self.model.results = solve_ivp(
            self.ode_system,
            [0, self.model.tmax],
            ic,
            t_eval=self.model.timeData,
            method='LSODA'
        )
        if not self.model.results.success:
            raise RuntimeError(f"ODE solver failed: {self.model.results.message}")

    def _postprocess_results(self):
        """Calculate derived quantities from ODE solution"""
        # Store positions
        self.model.bodyPosData = self.model.results.y[0]
        self.model.wheelPosData = self.model.results.y[2]

        # Road profile
        self.model.roadPosData = np.array([
            self.model.ymag if t > self.model.tramp
            else t * self.model.ymag / self.model.tramp
            for t in self.model.timeData
        ])

        # Calculate forces and accelerations
        self._calculate_forces()
        self._calculate_acceleration()

    def _calculate_forces(self):
        """Compute time history of all forces"""
        self.model.springForceData = []
        self.model.dashpotForceData = []
        self.model.tireForceData = []

        for i in range(len(self.model.timeData)):
            x1 = self.model.results.y[0][i]
            x2 = self.model.results.y[2][i]
            x1_dot = self.model.results.y[1][i]
            x2_dot = self.model.results.y[3][i]
            t = self.model.timeData[i]

            y_road = self.model.ymag if t > self.model.tramp else t * self.model.ymag / self.model.tramp

            self.model.springForceData.append(self.model.k1 * (x2 - x1))
            self.model.dashpotForceData.append(self.model.c1 * (x2_dot - x1_dot))
            self.model.tireForceData.append(self.model.k2 * (y_road - x2))

    def _calculate_acceleration(self):
        """Calculate body acceleration using finite differences"""
        vel = self.model.results.y[1]
        self.model.accelBodyData = np.zeros_like(vel)

        for i in range(len(vel)):
            if i == len(vel) - 1:
                h = self.model.timeData[i] - self.model.timeData[i - 1]
                self.model.accelBodyData[i] = (vel[i] - vel[i - 1]) / (9.81 * h)
            else:
                h = self.model.timeData[i + 1] - self.model.timeData[i]
                self.model.accelBodyData[i] = (vel[i + 1] - vel[i]) / (9.81 * h)

        self.model.accelMax = np.max(np.abs(self.model.accelBodyData))

    def OptimizeSuspension(self):
        """Optimization routine triggered by Optimize button"""
        try:
            self._load_parameters()
            x0 = np.array([self.model.k1, self.model.c1, self.model.k2])
            result = minimize(self.SSE, x0, method='Nelder-Mead')

            if result.success:
                self.model.k1, self.model.c1, self.model.k2 = result.x
                self.calculate()
                qtw.QMessageBox.information(None, "Optimization Complete",
                                            "New parameters:\n"
                                            f"k1: {self.model.k1:.1f} N/m\n"
                                            f"c1: {self.model.c1:.1f} N·s/m\n"
                                            f"k2: {self.model.k2:.1f} N/m")
            else:
                raise RuntimeError("Optimization failed: " + result.message)
        except Exception as e:
            qtw.QMessageBox.critical(None, "Optimization Error", str(e))

    def SSE(self, vals):
        """Sum of Squared Errors objective function for optimization"""
        k1, c1, k2 = vals
        self.model.k1 = k1
        self.model.c1 = c1
        self.model.k2 = k2

        try:
            self._solve_ode()
            self._postprocess_results()
        except:
            return float('inf')

        SSE = 0
        for i, t in enumerate(self.model.timeData):
            target = self.model.roadPosData[i]
            actual = self.model.bodyPosData[i]
            SSE += (actual - target) ** 2

        if self.chk_IncludeAccel and self.model.accelMax > self.model.accelLim:
            SSE += 10 * (self.model.accelMax - self.model.accelLim) ** 2

        return SSE

    def setupEventFilter(self, window):
        self.view.gv_Schematic.viewport().installEventFilter(window)

    def getPoints(self, t):
        """Get interpolated values at specific time t"""
        if self.model.results:
            idx = np.abs(self.model.timeData - t).argmin()
            return (
                self.model.results.y[0][idx],
                self.model.results.y[2][idx],
                self.model.roadPosData[idx],
                self.model.accelBodyData[idx]
            )
        return (0, 0, 0, 0)

    def doPlot(self):
        self.view.doPlot(self.model)

    def updateSchematic(self):
        self.view.updateSchematic()

    def setWidgets(self, w):
        self.chk_IncludeAccel = w

    def getZoom(self):
        """Get current zoom level from view"""
        return self.view.getZoom()

    def setZoom(self, zoom):
        """Update zoom level in view"""
        self.view.setZoom(zoom)



if __name__ == '__main__':
    # Test code
    model = CarModel()
    controller = CarController((None, None))
    controller.calculate()

