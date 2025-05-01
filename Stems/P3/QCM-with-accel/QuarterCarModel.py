# region imports
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import numpy as np
import math
from PyQt5 import QtWidgets as qtw, QtCore as qtc, QtGui as qtg
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure


class CarModel:
    def __init__(self):
        self.results = None
        self.timeData = np.linspace(0, 3, 2000)
        self.springForceData = []
        self.dashpotForceData = []
        self.tireForceData = []
        self.bodyPosData = []
        self.wheelPosData = []
        self.roadPosData = []
        self.accelBodyData = []
        self.tmax = 3.0
        self.tramp = 1.0
        self.angrad = 0.1
        self.ymag = 0.1515
        self.yangdeg = 45.0
        self.m1 = 450
        self.m2 = 20
        self.c1 = 4500
        self.k1 = 15000
        self.k2 = 90000
        self.v = 120.0
        self.accelMax = 0
        self.accelLim = 1.5
        self.SSE = 0.0


class CarView:
    def __init__(self, args):
        self.input_widgets, self.display_widgets = args
        (self.gv_Schematic, self.chk_LogX, self.chk_LogY,
         self.chk_LogAccel, self.chk_ShowAccel, self.lbl_MaxMinInfo,
         self.layout_pos, self.layout_forces) = self.display_widgets

        # Position plot setup
        self.figure_pos = Figure()
        self.canvas_pos = FigureCanvasQTAgg(self.figure_pos)
        self.ax_pos = self.figure_pos.add_subplot()
        self.layout_pos.addWidget(NavigationToolbar2QT(self.canvas_pos))
        self.layout_pos.addWidget(self.canvas_pos)

        # Force plot setup
        self.figure_force = Figure()
        self.canvas_force = FigureCanvasQTAgg(self.figure_force)
        self.ax_force = self.figure_force.add_subplot()
        self.layout_forces.addWidget(NavigationToolbar2QT(self.canvas_force))
        self.layout_forces.addWidget(self.canvas_force)

        self.buildScene()

    def setupEventFilter(self, window):
        """Add event filter to the schematic view"""
        self.gv_Schematic.setMouseTracking(True)
        self.gv_Schematic.scene().installEventFilter(window)

    def setupCanvasMoveEvent(self, window):
        """Connect canvas motion events"""
        self.canvas_pos.mpl_connect("motion_notify_event", window.mouseMoveEvent_Canvas)

    def buildScene(self):
        self.scene = qtw.QGraphicsScene()
        self.scene.setSceneRect(-200, -200, 400, 400)
        self.gv_Schematic.setScene(self.scene)
        # ... (original schematic drawing code)

    def updateView(self, model=None):
        # ... (original view update code)
        self.doPlot(model)
        self.doForcePlot(model)

    def doPlot(self, model=None):
        self.ax_pos.clear()
        if model and model.results:
            t = model.timeData
            self.ax_pos.plot(t, model.bodyPosData, 'b-', label='Body')
            self.ax_pos.plot(t, model.wheelPosData, 'r-', label='Wheel')
            self.ax_pos.plot(t, model.roadPosData, 'k-', label='Road')
            # ... (original plotting code)
        self.canvas_pos.draw()

    def doForcePlot(self, model=None):
        self.ax_force.clear()
        if model and model.springForceData:
            t = model.timeData
            self.ax_force.plot(t, model.springForceData, 'b-', label='Spring')
            self.ax_force.plot(t, model.dashpotForceData, 'r-', label='Dashpot')
            self.ax_force.plot(t, model.tireForceData, 'g-', label='Tire')
            self.ax_force.set_xlabel('Time (s)')
            self.ax_force.set_ylabel('Force (N)')
            self.ax_force.legend()
            self.ax_force.grid(True)
        self.canvas_force.draw()


class CarController:
    def __init__(self, args):
        self.input_widgets, self.display_widgets = args
        self.model = CarModel()
        self.view = CarView(args)

    def calculateForces(self):
        self.model.springForceData = []
        self.model.dashpotForceData = []
        self.model.tireForceData = []

        k1 = self.model.k1
        c1 = self.model.c1
        k2 = self.model.k2

        for i in range(len(self.model.timeData)):
            x1 = self.model.results.y[0][i]
            x2 = self.model.results.y[2][i]
            x1dot = self.model.results.y[1][i]
            x2dot = self.model.results.y[3][i]

            # Road position
            t = self.model.timeData[i]
            y = self.model.ymag if t > self.model.tramp else t * self.model.ymag / self.model.tramp

            # Force calculations
            spring_force = k1 * (x2 - x1)
            dashpot_force = c1 * (x2dot - x1dot)
            tire_force = k2 * (y - x2)

            self.model.springForceData.append(spring_force)
            self.model.dashpotForceData.append(dashpot_force)
            self.model.tireForceData.append(tire_force)

    def calculate(self):
        # Original parameter loading
        self.model.m1 = float(self.input_widgets[0].text())
        # ... (parameter loading code)

        # Solve ODE
        v = 1000 * self.model.v / 3600
        self.model.angrad = math.radians(self.model.yangdeg)
        self.model.tramp = self.model.ymag / (math.sin(self.model.angrad) * v)

        ic = [0, 0, 0, 0]
        self.model.results = solve_ivp(self.ode_system, [0, self.model.tmax], ic,
                                       t_eval=self.model.timeData)

        # Calculate derived data
        self.model.bodyPosData = self.model.results.y[0]
        self.model.wheelPosData = self.model.results.y[2]
        self.model.roadPosData = [self.model.ymag if t > self.model.tramp
                                  else t * self.model.ymag / self.model.tramp
                                  for t in self.model.timeData]

        self.calculateForces()
        self.calcAccel()
        self.view.updateView(self.model)

    # ... (rest of original controller methods)

    def setWidgets(self, w):
        """
        Pass widgets to view for setup.
        :param w:
        :return:
        """
        self.view.setWidgets(w)
        self.chk_IncludeAccel=self.view.chk_IncludeAccel

    def setupCanvasMoveEvent(self, window):
        """
        Pass through to view.
        :param window:
        :return:
        """
        self.view.setupCanvasMoveEvent(window)

    def setupEventFilter(self, window):
        self.view.setupEventFilter(window=window)

    def setupCanvasMoveEvent(self, window):
        self.view.setupCanvasMoveEvent(window)

    def getZoom(self):
        """
        Pass request along to the view.
        :return:
        """
        return self.view.getZoom()

    def setZoom(self,val):
        """
        Pass request along to the view.
        :param val:
        :return:
        """
        self.view.setZoom(val=val)

    def updateSchematic(self):
        """
        Pass request along to the view.
        :return:
        """
        self.view.updateSchematic()

    def doCalc(self, doPlot=True, doAccel=True):
        """
        This solves the differential equations for the quarter car model.
        :param doPlot:
        :param doAccel:
        :return:
        """
        v = 1000 * self.model.v / 3600  # convert speed to m/s from kph
        self.model.angrad = self.model.yangdeg * math.pi / 180.0  # convert angle to radians
        self.model.tramp = self.model.ymag / (math.sin(self.model.angrad) * v)  # calculate time to traverse ramp

        self.model.timeData=np.logspace(np.log10(0.000001), np.log10(self.model.tmax), 2000)
        #self.model.timeData=np.linspace(0, self.model.tmax, 2000)
        self.model.roadPosData=[self.model.ymag if t>self.model.tramp else t*self.model.ymag/self.model.tramp for t in self.model.timeData]
        ic = [0, 0, 0, 0]
        # run odeint solver
        self.step=0
        self.model.results = solve_ivp(self.ode_system, t_span=[0,self.model.tmax], y0=ic, t_eval=self.model.timeData)
        if doAccel:
            self.calcAccel()
        if doPlot:
            self.doPlot()
        self.model.bodyPosData = self.model.results.y[0]
        self.model.wheelPosData = self.model.results.y[2]
        pass

    def calcAccel(self):
        """
        Calculate the acceleration in the vertical direction using the forward difference formula.
        """
        N=len(self.model.timeData)
        self.model.accelBodyData=np.zeros(shape=N)
        vel=self.model.results.y[1]
        for i in range(N):
            if i==N-1:
                h = self.model.timeData[i] - self.model.timeData[i - 1]
                self.model.accelBodyData[i]= (vel[i] - vel[i - 1]) / (9.81 * h)  # backward difference of velocity
            else:
                h = self.model.timeData[i + 1] - self.model.timeData[i]
                self.model.accelBodyData[i] = (vel[i + 1] - vel[i]) / (9.81 * h)  # forward difference of velocity
            # else:
            #     self.model.accel[i]=(vel[i+1]-vel[i-1])/(9.81*2.0*h)  # central difference of velocity
        self.model.accelMax=self.model.accelBodyData.max()
        return True

    def OptimizeSuspension(self):
        """
        Step 1:  set parameters based on GUI inputs by calling self.set(doCalc=False)
        Step 2:  make an initial guess for k1, c1, k2
        Step 3:  optimize the suspension
        :return:
        """
        #Step 1:
        #$JES MISSING CODE HERE$
        self.calculate(doCalc=False)
        #Step 2:
        #JES MISSING CODE HERE$
        x0=np.array([(self.model.mink1)*1.1, self.model.c1, (self.model.mink2)*1.1])
        #Step 3:
        #JES MISSING CODE HERE$
        answer=minimize(self.SSE,x0,method='Nelder-Mead')
        self.view.updateView(self.model)

    def SSE(self, vals, optimizing=True):
        """
        Calculates the sum of square errors between the contour of the road and the car body.
        :param vals:
        :param optimizing:
        :return:
        """
        k1, c1, k2=vals  #unpack the new values for k1, c1, k2
        self.model.k1=k1
        self.model.c1=c1
        self.model.k2=k2
        self.doCalc(doPlot=False)  #solve the odesystem with the new values of k1, c1, k2
        SSE=0
        for i in range(len(self.model.results.y[0])):
            t=self.model.timeData[i]
            y=self.model.results.y[0][i]
            if t < self.model.tramp:
                ytarget = self.model.ymag * (t / self.model.tramp)
            else:
                ytarget = self.model.ymag
            SSE+=(y-ytarget)**2

        #some penalty functions if the constants are too small
        if optimizing:
            if k1<self.model.mink1 or k1>self.model.maxk1:
                SSE+=100
            if c1<10:
                SSE+=100
            if k2<self.model.mink2 or k2>self.model.maxk2:
                SSE+=100
            o_IncludeAccel = self.chk_IncludeAccel.isChecked()
            # I'm overlaying a gradient in the acceleration limit that scales with distance from a target squared.
            if self.model.accelMax > self.model.accelLim and o_IncludeAccel:
                # need to soften suspension
                SSE+=10+10*(self.model.accelMax-self.model.accelLim)**2
        self.model.SSE=SSE
        return SSE

    def doPlot(self):
        self.view.doPlot(self.model)

    def animate(self, t):
        self.view.animate(self.model,t)

    def getPoints(self, t):
        return self.view.getPoints(self.model, t)
#endregion
#endregion

def main():
    QCM = CarController()
    QCM.doCalc()

if __name__ == '__main__':
    main()

