import numpy as np
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from configuration import *
from dsfit import *
import matplotlib.pyplot as plt
import time
from IPython import display
import importlib, sys

import matplotlib.pyplot as plt
qmm = QuantumMachinesManager()

class OPXexp():
    def __init__(self):

        super().__init__()
    
    def resspec(self,*args):
        fmin=5e6
        fmax=15e6
        df=.1e6
        fvec=np.arange(fmin,fmax,df)
        fr=174.2e6

        with program() as prog:
            f = declare(int)
            n = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_st = declare_stream()
            Q_st = declare_stream()

            with for_(n,0,n<1000,n+1):
                with for_(f,fmin,f<fmax,f+df):
                    update_frequency("resonator",f)
                    update_frequency("qubit",fr)
                    play("pi", "qubit")
                    align("qubit","resonator")
                    reset_phase('resonator')
                    measure("long_readout"*amp(0.2),'resonator',None,
                            demod.full('long_integW_cos',I,"out1"),
                            demod.full('long_integW_sin',Q,"out1")
                           )
                    save(I,I_st)
                    save(Q,Q_st)
                    wait(3000,'resonator')
            with stream_processing():
                I_st.buffer(len(fvec)).average().save("I")
                Q_st.buffer(len(fvec)).average().save("Q")

        return prog,fvec
    
    def plotIQ(self,I,Q,fvec):
        fig = plt.figure(figsize = (10,5))
        ax = fig.add_subplot(111)
        
        ax.plot(fvec,I,label='I')
        ax.plot(fvec,Q,label='Q')
        mags = np.sqrt(I**2+Q**2)

        ax.plot(fvec,mags,'g.')

        ax.axvline(fvec[np.argmax(mags)])
        p = fitlor(fvec,mags**2,showfit=False)
        ax.plot(fvec,np.sqrt(lorfunc(p,fvec)),'r-')

        nu_res = p[2] + 7.33e9
        Qfac = nu_res/p[3]/2

        return (fig.show(),print('fr=',float('%.5g' % nu_res),'Q=',float('%.4g' % Qfac)))
    
    def qubitspec2d(self,*args):
        fmin=172e6
        fmax=177e6
        df=0.02e6
        reset_time=15e3
        fvec=np.arange(fmin,fmax,df)

        amin=0.0
        amax=0.05
        da=0.005
        avec=np.arange(amin,amax,da)
        avgs=1000
        
        with program() as prog:
            fr = declare(int)
            n = declare(int)
            I = declare(fixed)
            a = declare(fixed)
            Q = declare(fixed)
            I_stream = declare_stream()
            Q_stream = declare_stream()

            with for_(n,0,n<avgs,n+1):
                with for_(a,amin,a<amax-da/2,a+da):
                    with for_(fr,fmin,fr<fmax,fr+df):
                        update_frequency("resonator",10.4e6)
                        update_frequency("qubit",fr)
                        wait(int(reset_time//4), "qubit")# wait for the qubit to relax, several T1s
                        play("saturation"*amp(a), "qubit",duration = 15e3)
                        align("qubit", "resonator")
                        #play("saturation"*amp(a), "qubit",duration = 5e3)
                        reset_phase('resonator')
                        measure("readout"*amp(1),'resonator',None,
                                demod.full('integW1',I,"out1"),
                                demod.full('integW2',Q,"out1")
                               )
                        save(I,I_stream)
                        save(Q,Q_stream)

            with stream_processing():
                I_stream.buffer(len(fvec)).buffer(len(avec)).average().save("I")
                Q_stream.buffer(len(fvec)).buffer(len(avec)).average().save("Q")
        
        return prog,fvec,avec
    
    def rabi2d(self,*args):
        dt = 1
        T_min = 1
        T_max = 2000//4
        times = np.arange(T_min, T_max, dt)*4

        fmin = 174e6-5e6
        fmax = 174e6+5e6
        df = 0.2e6
        fvec=np.arange(fmin,fmax,df)
        reset_time = 5*T1

        avgs = 1000
        with program() as prog:

            ##############################
            # declare real-time variables:
            ##############################

            n = declare(int)      # Averaging
            i = declare(int)      # Amplitudes
            t = declare(int)     #array of time delays
            f = declare(int)
            I = declare(fixed)
            Q = declare(fixed)

            I_st = declare_stream()
            Q_st = declare_stream()


            ###############
            # the sequence:
            ###############
            update_frequency("resonator",9.5e6)
            with for_(n, 0, n < avgs, n + 1):
                with for_(f, fmin,f<fmax, f+df):
                    update_frequency("qubit", f)
                    with for_(t, T_min, t < T_max, t + dt):
                        wait(int(reset_time//4), "qubit")
                        play("saturation", "qubit",duration = t)
                        align("qubit", "resonator")
                        reset_phase('resonator')
                        measure("readout"*amp(1),'resonator',None,
                                demod.full('integW1',I,"out1"),
                                demod.full('integW2',Q,"out1")
                               )
                        save(I,I_st)
                        save(Q,Q_st)

#                        """Play a ge pi pulse and then readout"""
#                 wait(reset_time // 4, "qubit")
#                 play("pi", "qubit")
#                 align('qubit','resonator')
#                 measure("readout"*amp(0.2),'resonator',None,
#                         demod.full('integW1',Ie,"out1"),
#                         demod.full('integW2',Qe,"out1")
#                                )

#                 save(Ie,Ie_st)
#                 save(Qe,Qe_st)

#                         #  """Just readout without playing anything"""
#                 wait(reset_time // 4, "qubit")
#                 align('qubit','resonator')
#                 measure("readout"*amp(0.2),'resonator',None,
#                         demod.full('integW1',Ig,"out1"),
#                         demod.full('integW2',Qg,"out1"))
#                 save(Ig,Ig_st)
#                 save(Qg,Qg_st)

#                align("qubit", "resonator")
            with stream_processing():
                I_st.buffer(len(times)).buffer(len(fvec)).average().save('I')
                Q_st.buffer(len(times)).buffer(len(fvec)).average().save('Q')

        return prog,fvec,times

    def plot2d(self,I,Q,x1,x2):
        fig = plt.figure(figsize = (10,5))
        ax = fig.add_subplot(111)
        
        mags = np.sqrt(I**2+Q**2)

        ax.pcolormesh(x1,x2,mags)
        return ax
    
    def T1(self,*args):
        dt = 20
        T_min = 1
        T_max = 10000
        times = np.arange(T_min, T_max, dt)*4
        T1=10.5e3
        reset_time = 5*T1

        avgs = 50000
        with program() as prog:

            ##############################
            # declare real-time variables:
            ##############################

            n = declare(int)      # Averaging
            i = declare(int)      # Amplitudes
            t = declare(int) #array of time delays
            f = declare(int)
            I = declare(fixed)
            Q = declare(fixed)

            I_st = declare_stream()
            Q_st = declare_stream()

            ###############
            # the sequence:
            ###############
            update_frequency("resonator",10.2e6)
            update_frequency("qubit", 174.2e6)
            with for_(n, 0, n < avgs, n + 1):
                with for_(t, T_min, t < T_max, t + dt):
                    wait(int(reset_time//4), "qubit")
                    play("saturation", "qubit",duration=193//4) # pi pulse with saturation
                    wait(t)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(1),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))
                    save(I,I_st)
                    save(Q,Q_st)

                wait(int(reset_time//4), "qubit")
                play("saturation", "qubit",duration=193//4) # pi pulse with saturation
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(1),'resonator',None, #amp=0.45, f=10.2 looks best for low power (resonator punched out) readout
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1")
                               )
                save(I,I_st)
                save(Q,Q_st)

                wait(int(reset_time//4), "qubit")
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(1),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))

                save(I,I_st)
                save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)+2).average().save('I')
                Q_st.buffer(len(times)+2).average().save('Q')

        return prog
    
    def ramsey2d(self, *args):
        dt = 1
        Tmin = 0
        Tmax = 1000
        times = np.arange(Tmin, Tmax, dt)*4

        fr=174.2e6

        pmin = 0
        pmax = 1
        dp = 0.1
        pvec = np.arange(pmin, pmax, dp)
        reset_time = 5*T1

        avgs = 1000
        with program() as prog:

            ##############################
            # declare real-time variables:
            ##############################

            n = declare(int)      # Averaging
            i = declare(int)      # Amplitudes
            t = declare(int) #array of time delays
            p = declare(fixed)
            I = declare(fixed)
            Q = declare(fixed)

            I_st = declare_stream()
            Q_st = declare_stream()


            ###############
            # the sequence:
            ###############
            update_frequency("qubit", 174.2e6)
            update_frequency("resonator", 10.2e6)

            with for_(n, 0, n < avgs, n + 1):

                with for_(p, pmin, p < pmax, p + dp):
                    assign(phi,0)
                    with for_(t, 0, t < T_max, t + dt):
                        wait(int(reset_time//4), "qubit")
                       # play("pi2", "qubit")
                        play("saturation","qubit",duration=193//8)
                        wait(t)
                        frame_rotation_2pi(phi,'qubit')
                      #  play("pi2", "qubit")
                        play("saturation","qubit",duration=193//8)
                        align("qubit", "resonator")
                        reset_phase('resonator')
                        measure("readout"*amp(1),'resonator',None,
                                demod.full('integW1',I,"out1"),
                                demod.full('integW2',Q,"out1")
                               )
                        assign(phi,phi+p*0.01*dt/4)
                        save(I,I_st)
                        save(Q,Q_st)
                        reset_frame('qubit')
                                #  """Just readout without playing anything"""
                    wait(int(reset_time // 4), "qubit")

                    play("saturation","qubit",duration=193//4)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(1),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1")
                                   )
                    save(I,I_st)
                    save(Q,Q_st)

                    wait(int(reset_time//4), "qubit")
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(1),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1"))

                    save(I,I_st)
                    save(Q,Q_st)
            with stream_processing():
                I_st.buffer(len(times)+2).buffer(len(pvec)).average().save('I')
                Q_st.buffer(len(times)+2).buffer(len(pvec)).average().save('Q')
            return prog
        
    def ramsey(self, *args):
        dt = 2
        T_min = 1
        T_max = 2000//4
        times = np.arange(T_min, T_max, dt)*4

        reset_time = 5*T1

        dphi_min = -0.05
        dphi_max = 0.05
        ddphi = 0.001
        dphi_vec = np.arange(dphi_min, dphi_max + ddphi/2, ddphi)

        avgs = 5000
        with program() as prog:

            ##############################
            # declare real-time variables:
            ##############################

            n = declare(int)      # Averaging
            i = declare(int)      # Amplitudes
            t = declare(int) #array of time delays
            f = declare(int)
            phi = declare(fixed)
            dphi = declare(fixed)
            I = declare(fixed)
            Q = declare(fixed)

            I_st = declare_stream()
            Q_st = declare_stream()

            update_frequency("resonator",10.2e6)
            update_frequency("qubit", 174.2e6)

            with for_(n, 0, n < avgs, n + 1):
             #   with for_(dphi, dphi_min, dphi < dphi_max + ddphi/2, dphi + ddphi):
                assign(phi,0)
                with for_(t, T_min, t < T_max, t + dt):
                    wait(int(reset_time//4), "qubit")
                    play("saturation","qubit",duration=193//8) #pi/2
                    wait(t)
                    frame_rotation_2pi(phi,'qubit')
                    play("saturation","qubit",duration=193//8)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(0.45),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1"))
                    assign(phi,phi+0.01*dt)
                    save(I,I_st)
                    save(Q,Q_st)
                    reset_frame('qubit')

                wait(int(reset_time // 4), "qubit")
                play("saturation","qubit",duration=193//4)
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.45),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1")
                               )
                save(I,I_st)
                save(Q,Q_st)

                wait(int(reset_time //4), "qubit")
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.45),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))

                save(I,I_st)
                save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)+2).average().save('I')
                Q_st.buffer(len(times)+2).average().save('Q')
        return prog

    def spin_echo(self,phase=0,*args):    
        dt = 5
        T_min = 2
        T_max = 4000//4
        times = np.arange(T_min, T_max, dt)*4
        T1=9e3
        reset_time = 5*T1

        avgs = 50000
        with program() as prog:

            ##############################
            # declare real-time variables:
            ##############################

            n = declare(int)      # Averaging
            i = declare(int)      # Amplitudes
            t = declare(int) #array of time delays
            f = declare(int)
            phi = declare(fixed)
            dphi = declare(fixed)
            I = declare(fixed)
            Q = declare(fixed)

            I_st = declare_stream()
            Q_st = declare_stream()

            ###############
            # the sequence:
            ###############
            update_frequency("resonator",10.2e6)
            update_frequency("qubit", 174.2e6)

            with for_(n, 0, n < avgs, n + 1):

                with for_(t, T_min, t < T_max, t + dt):
                    wait(int(reset_time//4), "qubit")
                    play("saturation", "qubit",duration=193//8) #pi/2
                    wait(t/2)
                    frame_rotation_2pi(phase,'qubit')
                    play("saturation", "qubit",duration=193//4) #pi
                    reset_frame('qubit')
                    wait(t/2)
                    frame_rotation_2pi(phi,'qubit')
                    play("saturation", "qubit",duration=193//8) #pi/2 with phase rotation
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(1),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1"))
                    assign(phi,phi+0.01*dt)
                    save(I,I_st)
                    save(Q,Q_st)
                    reset_frame('qubit')

                wait(int(reset_time // 4), "qubit")
                play("saturation", "qubit",duration=193//4)
                reset_phase('resonator')
                align("qubit", "resonator")
                measure("readout"*amp(1),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1")
                               )
                save(I,I_st)
                save(Q,Q_st)

                wait(int(reset_time //4), "qubit")
                reset_phase('resonator')
                align("qubit", "resonator")
                measure("readout"*amp(1),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))

                save(I,I_st)
                save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)+2).average().save('I')
                Q_st.buffer(len(times)+2).average().save('Q')
            
            return prog
        
    def ef_spec(self,*args):
        fmin=10e6
        fmax=30e6
        df=0.1e6
        reset_time=15e3
        fvec=np.arange(fmin,fmax,df)

        avgs=50e3
        with program() as prog:
            fr = declare(int)
            n = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_stream = declare_stream()
            Q_stream = declare_stream()

            with for_(n,0,n<avgs,n+1):
                with for_(fr,fmin,fr<fmax,fr+df):
                    update_frequency("resonator",10.2e6)
                    update_frequency("qubit",174.2e6)
                    wait(int(reset_time//4), "qubit")# wait for the qubit to relax, several T1s
                    play("saturation","qubit",duration=193//4)
                    update_frequency("qubit",fr)
                    play("saturation", "qubit",duration = 1e3)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(1),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1")
                           )
                    save(I,I_stream)
                    save(Q,Q_stream)

            with stream_processing():
                I_stream.buffer(len(fvec)).average().save("I")
                Q_stream.buffer(len(fvec)).average().save("Q")
                
        return prog