from __future__ import division
from neuron import h
from neuron import gui
import matplotlib.pyplot as plt
import numpy as np
plt.ion()

# create model
h("create soma")
h.soma.L    = 10
h.soma.diam = 10
h.soma.Ra   = 100
h.soma.insert('pas')   # add passive properties  
h.soma.g_pas = 1/10000 # set the specific membrane resistance to 10000 ohm*cm^2

# add active conductances (the channels [mod files] are from Mainen and Sejnowski 1996)
h.soma.insert('kv') # add potassium channel
h.soma.gbar_kv = 2000 # set the potassium conductance

h.soma.insert('na') # add sodium channel
h.soma.gbar_na = 8000 # set the sodium conductance
h.celsius = 30


# current clamp
stim = h.IClamp(h.soma(0.5))
stim.amp   = 0.007446
stim.delay = 250
stim.dur   = 1000

# record voltage of some and injected current
# and the time
soma_v = h.Vector()
soma_v.record(h.soma(0.5)._ref_v)

stim_current = h.Vector()
stim_current.record(stim._ref_i)

t = h.Vector()
t.record(h._ref_t)


# run simulation
h.tstop = 1500 # set the simulation time
h.v_init = -70
h.run()

# plot the injected current and the voltage response
f, (ax0, ax1) = plt.subplots(2,1, figsize=(4,3), gridspec_kw = {'height_ratios':[3, 1]})
ax0.plot(t,soma_v, 'k')
ax1.plot(t,stim_current, 'gray')

ax0.set_ylabel('Voltage (mV)')
ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.get_xaxis().set_visible(False)


ax1.plot([50,50],[0.01,0.02],'k')
ax1.text(80,0.015,'10pA',va='center')
ax1.set_ylabel('I (nA)')
ax1.set_xlabel('t (ms)')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.get_yaxis().set_visible(False)
plt.tight_layout()
