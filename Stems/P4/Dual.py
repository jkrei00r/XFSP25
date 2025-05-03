#region imports
from Air import *
from matplotlib import pyplot as plt
from PyQt5 import QtWidgets as qtw
import numpy as np
import sys
#endregion

#region class definitions
class dualCycleModel:
    def __init__(self, p_initial=100000, v_cylinder=0.003, t_initial=300, pressure_ratio=1.5, cutoff_ratio=1.2, compression_ratio=18, name='Air Standard Dual Cycle'):
        self.units = units()
        self.air = air()
        self.p_initial = p_initial
        self.t_initial = t_initial
        self.compression_ratio = compression_ratio
        self.pressure_ratio = pressure_ratio  # P3/P2
        self.cutoff_ratio = cutoff_ratio      # V4/V3
        self.v_cylinder = v_cylinder
        self.upperCurve = StateDataForPlotting()
        self.lowerCurve = StateDataForPlotting()
        self.calculated = False
        self.cycleType = 'dual'
        self.set()

    def set(self, SI=True, T_0=300, P_0=100000, V_0=0.003, PR=1.5, CR=1.2, r=18):
        self.units.set(SI=SI)
        self.t_initial = T_0 if SI else self.units.T_RtoK(T_0)
        self.p_initial = P_0 if SI else P_0 / self.units.CF_P
        self.pressure_ratio = PR
        self.cutoff_ratio = CR
        self.compression_ratio = r
        self.v_cylinder = V_0 if SI else V_0 / self.units.CF_V

        # State 1: Initial state
        self.State1 = self.air.set(P=self.p_initial, T=self.t_initial)
        # State 2: Isentropic compression
        self.State2 = self.air.set(v=self.State1.v / self.compression_ratio, s=self.State1.s)
        # State 3: Constant volume heat addition (P3 = P2 * pressure_ratio)
        self.State3 = self.air.set(P=self.State2.P * self.pressure_ratio, v=self.State2.v)
        # State 4: Constant pressure heat addition (V4 = V3 * cutoff_ratio)
        self.State4 = self.air.set(P=self.State3.P, v=self.State3.v * self.cutoff_ratio)
        # State 5: Isentropic expansion
        self.State5 = self.air.set(v=self.State1.v, s=self.State4.s)

        # Calculate work and heat
        self.W_compression = self.air.n * (self.State2.u - self.State1.u)
        self.W_power = self.air.n * ((self.State4.u - self.State5.u) + self.State3.P * (self.State4.v - self.State3.v))
        self.Q_in = self.air.n * ((self.State3.u - self.State2.u) + (self.State4.h - self.State3.h))
        self.Q_out = self.air.n * (self.State5.u - self.State1.u)
        self.efficiency = 100 * (self.W_power - self.W_compression) / self.Q_in if self.Q_in != 0 else 0

        self.calculated = True
        self.buildDataForPlotting()

    def buildDataForPlotting(self):
        self.upperCurve.clear()
        self.lowerCurve.clear()
        a = air()

        # Process 2-3: Constant volume
        for P in np.linspace(self.State2.P, self.State3.P, 30):
            state = a.set(P=P, v=self.State2.v)
            self.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

        # Process 3-4: Constant pressure
        for v in np.linspace(self.State3.v, self.State4.v, 30):
            state = a.set(P=self.State3.P, v=v)
            self.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

        # Process 4-5: Isentropic expansion
        for v in np.linspace(self.State4.v, self.State5.v, 30):
            state = a.set(v=v, s=self.State4.s)
            self.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

        # Process 5-1: Constant volume heat rejection
        for T in np.linspace(self.State5.T, self.State1.T, 30):
            state = a.set(T=T, v=self.State5.v)
            self.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

        # Process 1-2: Isentropic compression
        for v in np.linspace(self.State1.v, self.State2.v, 30):
            state = a.set(v=v, s=self.State1.s)
            self.lowerCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

class dualCycleController:
    def __init__(self, model=None, ax=None):
        self.model = model if model else dualCycleModel()
        self.view = dualCycleView()
        self.view.ax = ax

    def calc(self):
        try:
            T0 = float(self.view.le_TLow.text())
            P0 = float(self.view.le_P0.text())
            V0 = float(self.view.le_V0.text())
            PR = float(self.view.le_THigh.text())  # Pressure ratio
            CR = float(self.view.le_CR.text())     # Cutoff ratio
            r = 18  # Default compression ratio
            self.set(T_0=T0, P_0=P0, V_0=V0, PR=PR, CR=CR, r=r, SI=self.view.rdo_Metric.isChecked())
        except ValueError:
            pass

    def set(self, T_0=300, P_0=100000, V_0=0.003, PR=1.5, CR=1.2, r=18, SI=True):
        self.model.set(SI=SI, T_0=T_0, P_0=P_0, V_0=V_0, PR=PR, CR=CR, r=r)
        self.updateView()

    def updateView(self):
        self.view.updateView(self.model)

    def setWidgets(self, w):
        [self.view.lbl_THigh, self.view.lbl_TLow, self.view.lbl_P0, self.view.lbl_V0, self.view.lbl_CR,
         self.view.le_THigh, self.view.le_TLow, self.view.le_P0, self.view.le_V0, self.view.le_CR,
         self.view.le_T1, self.view.le_T2, self.view.le_T3, self.view.le_T4,
         self.view.lbl_T1Units, self.view.lbl_T2Units, self.view.lbl_T3Units, self.view.lbl_T4Units,
         self.view.le_PowerStroke, self.view.le_CompressionStroke, self.view.le_HeatAdded, self.view.le_Efficiency,
         self.view.lbl_PowerStrokeUnits, self.view.lbl_CompressionStrokeUnits, self.view.lbl_HeatInUnits,
         self.view.rdo_Metric, self.view.cmb_Abcissa, self.view.cmb_Ordinate,
         self.view.chk_LogAbcissa, self.view.chk_LogOrdinate, self.view.ax, self.view.canvas] = w

class dualCycleView:
    def __init__(self):
        self.lbl_THigh = qtw.QLabel()
        self.le_THigh = qtw.QLineEdit()
        self.lbl_CR = qtw.QLabel()
        self.le_CR = qtw.QLineEdit()
        self.ax = None
        self.canvas = None

    def updateView(self, cycle):
        if not cycle.calculated:
            return

        # Update labels
        self.lbl_THigh.setText('Pressure Ratio (P3/P2)')
        self.lbl_CR.setText('Cutoff Ratio (rc)')

        # Update plot
        self.plot_cycle_XY(cycle)

    def plot_cycle_XY(self, cycle):
        self.ax.clear()
        # Plotting logic similar to Otto/Diesel
        # (Implement based on cycle.upperCurve and cycle.lowerCurve)
        self.canvas.draw()
#endregion
