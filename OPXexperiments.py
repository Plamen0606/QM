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
        fmin=55e6
        fmax=65e6
        df=0.1e6
        fvec=np.arange(fmin,fmax,df)
        fr=174.2e6

        with program() as prog:
            f = declare(int)
            n = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_st = declare_stream()
            Q_st = declare_stream()

            with for_(n,0,n<5000,n+1):
                with for_(f,fmin,f<fmax,f+df):
                    update_frequency("resonator",f)
                    align("qubit","resonator")
                  #  reset_phase('resonator')
                    measure("readout"*amp(1),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1")
                           )
                    save(I,I_st)
                    save(Q,Q_st)
                    wait(3000,'resonator')
            with stream_processing():
                I_st.buffer(len(fvec)).average().save("I")
                Q_st.buffer(len(fvec)).average().save("Q")

        return prog,fvec
    def resspec_amp(self,*args):
        fmin = 55e6
        fmax = 65e6
        df = .2e6
        fvec = np.arange(fmin, fmax, df)

        amin = 0.01
        amax = 1
        da = 0.1
        avec = np.arange(amin, amax+da/2, da)
        reset_time=15e3

        with program() as prog:
            fr = declare(int)
            n = declare(int)
            a = declare(fixed)
            I1 = declare(fixed)
            Q1 = declare(fixed)
            I1_stream = declare_stream()
            Q1_stream = declare_stream()

            with for_(n, 0, n < 2000, n + 1):
                with for_(fr, fmin, fr < fmax, fr + df):
                    with for_(a, amin, a < amax + da / 2, a + da):
                        update_frequency("resonator", fr)
                        # align("Vsource", "resonator")
                        #  play("saturation"*amp(a), "resonator",duration = 20e3) #0.8 V for amp=1
                        #reset_phase('resonator')
                        measure("readout" * amp(a), 'resonator', None,
                                demod.full('integW1', I1, "out1"),
                                demod.full('integW2', Q1, "out1")
                                )
                        save(I1, I1_stream)
                        save(Q1, Q1_stream)
                        wait(int(reset_time//4), 'resonator')
            with stream_processing():
                I1_stream.buffer(len(avec)).buffer(len(fvec)).average().save("I")
                Q1_stream.buffer(len(avec)).buffer(len(fvec)).average().save("Q")
        return prog, fvec, avec

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

        return (fig.show(),print('fr=',float('%.5g' % nu_res),print('frIF=',float('%.5g' % p[2]),'Q=',float('%.4g' % Qfac))))

    def qubitspec2d(self,*args):
        fmin=-27e6
        fmax=-25.5e6
        df=0.1e6
        reset_time=10e3
        fvec=np.arange(fmin,fmax,df)

        amin=0.6
        amax=1
        da=0.05
        avec=np.arange(amin,amax+da/2,da)
        avgs=20000

        with program() as prog:
            fr = declare(int)
            n = declare(int)
            I = declare(fixed)
            a = declare(fixed)
            Q = declare(fixed)
            I_stream = declare_stream()
            Q_stream = declare_stream()

            with for_(n,0,n<avgs,n+1):
                with for_(a,amin,a<amax+da/2,a+da):
                    with for_(fr,fmin,fr<fmax,fr+df):
                        update_frequency("resonator",59.8e6)
                        update_frequency("qubit",fr)
                        wait(int(reset_time//4), "qubit")# wait for the qubit to relax, several T1s
                        #play("saturation"*amp(1), "qubit",duration = 15e3)
                        align("qubit", "resonator")
                        play("saturation"*amp(a), "qubit",duration = 1e3//4)
                        reset_phase('resonator')
                        measure("readout"*amp(0.6),'resonator',None,
                                demod.full('integW1',I,"out1"),
                                demod.full('integW2',Q,"out1")
                               )
                        save(I,I_stream)
                        save(Q,Q_stream)

            with stream_processing():
                I_stream.buffer(len(fvec)).buffer(len(avec)).average().save("I")
                Q_stream.buffer(len(fvec)).buffer(len(avec)).average().save("Q")
        
        return prog,fvec+4.1e9,avec

    def qubitspec(self, *args):
        fmin = 400e6-27e6
        fmax = 400e6-25.5e6
        df = 0.1e6
        reset_time = 2e3
        fvec = np.arange(fmin, fmax, df)

        avgs = 10000

        with program() as prog:
            fr = declare(int)
            n = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_stream = declare_stream()
            Q_stream = declare_stream()

            with for_(n, 0, n < avgs, n + 1):
                with for_(fr, fmin, fr < fmax, fr + df):
                    update_frequency("resonator", 59.8e6)
                    update_frequency("qubit", fr)
                    wait(int(reset_time // 4), "qubit")  # wait for the qubit to relax, several T1s
                   # play("saturation" * amp(1), "qubit", duration=5e4, chirp=(1e6,'MHz/sec'))
                    play("saturation"*amp(1), "qubit",duration = 1e3//4)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout" * amp(0.6), 'resonator', None,
                            demod.full('integW1', I, "out1"),
                            demod.full('integW2', Q, "out1")
                            )
                    save(I, I_stream)
                    save(Q, Q_stream)

            with stream_processing():
                I_stream.buffer(len(fvec)).average().save("I")
                Q_stream.buffer(len(fvec)).average().save("Q")

        return prog, fvec+3.7e9

    def qubitspec_gauss(self, *args):
        fmin = 400e6 - 27e6
        fmax = 400e6 - 25.5e6
        df = 0.1e6
        reset_time = 2e3
        fvec = np.arange(fmin, fmax, df)

        avgs = 10000

        with program() as prog:
            fr = declare(int)
            n = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_stream = declare_stream()
            Q_stream = declare_stream()

            with for_(n, 0, n < avgs, n + 1):
                with for_(fr, fmin, fr < fmax, fr + df):
                    update_frequency("resonator", 59.8e6)
                    update_frequency("qubit", fr)
                    wait(int(reset_time // 4), "qubit")  # wait for the qubit to relax, several T1s
                    # play("saturation" * amp(1), "qubit", duration=5e4, chirp=(1e6,'MHz/sec'))
                    play("gaussian" * amp(1), "qubit", len=1e3 // 4)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout" * amp(0.6), 'resonator', None,
                            demod.full('integW1', I, "out1"),
                            demod.full('integW2', Q, "out1")
                            )
                    save(I, I_stream)
                    save(Q, Q_stream)

            with stream_processing():
                I_stream.buffer(len(fvec)).average().save("I")
                Q_stream.buffer(len(fvec)).average().save("Q")

        return prog, fvec + 3.7e9

    def qubitspec_V(self, *args):
        fmin=-27e6
        fmax=-25.5e6
        df=0.02e6
        reset_time=10e3

        fvec = np.arange(fmin, fmax, df)

        avgs = 100000

        amin=0
        amax=1
        da=0.05
        avec=np.arange(amin,amax,da)

        with program() as prog:
            fr = declare(int)
            n = declare(int)
            a= declare(fixed)
            I = declare(fixed)
            Q = declare(fixed)
            I_stream = declare_stream()
            Q_stream = declare_stream()
            update_frequency("resonator", 59.8e6)

            with for_(n, 0, n < avgs, n + 1):
                with for_(a,amin,a<amax+da/2,a+da):
                    align("Vsource","resonator", "qubit")
                    play("CW" * amp(a), "Vsource", duration=int(len(fvec)*(reset_time+5e3)// 4))  # 0.8 V for amp=
                    with for_(fr, fmin, fr < fmax, fr + df):
                        update_frequency("qubit", fr)
                        wait(int(reset_time // 4), "qubit")  # wait for the qubit to relax, several T1s
                        play("saturation", "qubit",duration = 5e3//4)
                        align("resonator", "qubit")
                        reset_phase('resonator')
                        measure("readout" * amp(0.6), 'resonator', None,
                                demod.full('integW1', I, "out1"),
                                demod.full('integW2', Q, "out1")
                                )
                        save(I, I_stream)
                        save(Q, Q_stream)

            with stream_processing():
                I_stream.buffer(len(fvec)).buffer(len(avec)).average().save("I")
                Q_stream.buffer(len(fvec)).buffer(len(avec)).average().save("Q")

        return prog, fvec+4.1e9, avec
    
    def rabi2d(self,*args):
        dt = 2
        T_min = 1
        T_max = 300//4
        times = np.arange(T_min, T_max, dt)*4

        T1=2e3

        fmin = 400e6 - 27e6
        fmax = 400e6 - 25.5e6
        df = 0.02e6
       # reset_time = 2e3
        fvec=np.arange(fmin,fmax,df)
        reset_time = 5*T1

        avgs = 10000
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
            update_frequency("resonator",59.8e6)
            with for_(n, 0, n < avgs, n + 1):
                with for_(f, fmin,f<fmax, f+df):
                    update_frequency("qubit", f)
                    with for_(t, T_min, t < T_max, t + dt):
                        wait(int(reset_time//4), "qubit")
                        play("saturation", "qubit",duration = t)
                        align("qubit", "resonator")
                        reset_phase('resonator')
                        measure("readout"*amp(0.6),'resonator',None,
                                demod.full('integW1',I,"out1"),
                                demod.full('integW2',Q,"out1")
                               )
                        save(I,I_st)
                        save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)).buffer(len(fvec)).average().save('I')
                Q_st.buffer(len(times)).buffer(len(fvec)).average().save('Q')

        return prog,fvec+3.7e9,times

    def rabi1d(self, *args):
        dt = 2
        T_min = 1
        T_max = 600 // 4
        times = np.arange(T_min, T_max, dt) * 4

        T1 = 10e3
        reset_time = 5 * T1

        avgs = 20000
        with program() as prog1:
            ##############################
            # declare real-time variables:
            ##############################

            n = declare(int)  # Averaging
            t = declare(int)  # array of time delays
            I = declare(fixed)
            Q = declare(fixed)

            I_st = declare_stream()
            Q_st = declare_stream()

            ###############
            # the sequence:
            ###############
            update_frequency("resonator", 59.4e6)
            update_frequency("qubit", -28e6)

            with for_(n, 0, n < avgs, n + 1):
                with for_(t, T_min, t < T_max, t + dt):
                    wait(int(reset_time // 4), "qubit")
                    play("saturation", "qubit", duration=t)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout" * amp(0.6), 'resonator', None,
                            demod.full('integW1', I, "out1"),
                            demod.full('integW2', Q, "out1")
                            )
                    save(I, I_st)
                    save(Q, Q_st)

                """Play a ge pi pulse and then readout"""
                wait(int(reset_time // 4), "qubit")
                play("saturation", "qubit",duration=40//4)
                align('qubit','resonator')
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))
                save(I,I_st)
                save(Q,Q_st)

                """Just readout without playing anything"""
                wait(int(reset_time // 4), "qubit")
                align('qubit','resonator')
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))
                save(I,I_st)
                save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)+2).average().save('I')
                Q_st.buffer(len(times)+2).average().save('Q')

        return prog1, times

    def plot2d(self,I,Q,x1,x2):
        fig = plt.figure(figsize = (10,5))
        ax = fig.add_subplot(111)
        
        mags = np.sqrt(I**2+Q**2)

        ax.pcolormesh(x1,x2,mags)
        return ax
    
    def T1(self,*args):
        dt = 20//4
        T_min = 1
        T_max = 3000//4
        times = np.arange(T_min, T_max, dt)*4
        T1=10.5e3
        reset_time = 5*T1

        avgs = 1000
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
            update_frequency("resonator",59.8e6)
            update_frequency("qubit", -26.2e6)
            with for_(n, 0, n < avgs, n + 1):
                with for_(t, T_min, t < T_max, t + dt):
                    wait(int(reset_time//4), "qubit")
                    play("saturation", "qubit",duration=40//4) # pi pulse with saturation
                    wait(t)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))
                    save(I,I_st)
                    save(Q,Q_st)

                wait(int(reset_time//4), "qubit")
                play("saturation", "qubit", duration=40// 4)  # pi pulse with saturation
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None, #amp=0.45, f=10.2 looks best for low power (resonator punched out) readout
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1")
                               )
                save(I,I_st)
                save(Q,Q_st)

                wait(int(reset_time//4), "qubit")
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))

                save(I,I_st)
                save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)+2).average().save('I')
                Q_st.buffer(len(times)+2).average().save('Q')

        return prog,times
    
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
        T_max = 500//4
        times = np.arange(T_min, T_max, dt)*4
        T1=10e3
        reset_time = 1e6

        dphi_min = -0.05
        dphi_max = 0.05
        ddphi = 0.001
        dphi_vec = np.arange(dphi_min, dphi_max + ddphi/2, ddphi)

        avgs = 1000000
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

            update_frequency("resonator",59.8e6)
            update_frequency("qubit", -24.5e6)

            with for_(n, 0, n < avgs, n + 1):
             #   with for_(dphi, dphi_min, dphi < dphi_max + ddphi/2, dphi + ddphi):
                assign(phi,0)
                with for_(t, T_min, t < T_max, t + dt):
                    wait(int(reset_time//4), "qubit")
                    play("saturation","qubit",duration=40//8) #pi/2
                    wait(t)
                    frame_rotation_2pi(phi,'qubit')
                    play("saturation","qubit",duration=40//8)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(0.6),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1"))
                    assign(phi,phi+0.1*dt)
                    save(I,I_st)
                    save(Q,Q_st)
                    reset_frame('qubit')

                wait(int(reset_time // 4), "qubit")
                play("saturation","qubit",duration=40//4)
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1")
                               )
                save(I,I_st)
                save(Q,Q_st)

                wait(int(reset_time //4), "qubit")
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))

                save(I,I_st)
                save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)+2).average().save('I')
                Q_st.buffer(len(times)+2).average().save('Q')
        return prog,times

    def spin_echo(self,phase=0,*args):    
        dt = 2
        T_min = 2
        T_max = 500//4
        times = np.arange(T_min, T_max, dt)*4
        T1=10e3
        reset_time = 5*T1
        phase=0.25

        avgs = 100000
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
            update_frequency("resonator",59.8e6)
            update_frequency("qubit", -24.5e6)

            with for_(n, 0, n < avgs, n + 1):
                assign(phi,0)
                with for_(t, T_min, t < T_max, t + dt):
                    wait(int(reset_time//4), "qubit")

                    play("saturation", "qubit",duration=40//8) #pi/2
                    wait(t/2)

                    frame_rotation_2pi(phase,'qubit')
                    play("saturation", "qubit",duration=40//4) #pi
                    reset_frame('qubit')
                    wait(t/2)

                    frame_rotation_2pi(phi,'qubit')
                    play("saturation", "qubit",duration=40//8) #pi/2 with phase rotation

                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(0.6),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1"))
                    assign(phi,phi+0.1*dt)
                    save(I,I_st)
                    save(Q,Q_st)
                    reset_frame('qubit')

                wait(int(reset_time // 4), "qubit")
                play("saturation", "qubit",duration=40//4)
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1")
                               )
                save(I,I_st)
                save(Q,Q_st)

                wait(int(reset_time //4), "qubit")
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))

                save(I,I_st)
                save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)+2).average().save('I')
                Q_st.buffer(len(times)+2).average().save('Q')
            
            return prog,times
        
    def ef_spec(self,pi=True,*args):
        fmin=10e6
        fmax=300e6
        df=0.2e6
        reset_time=10e3
        fvec=np.arange(fmin,fmax,df)

        avgs=5e3
        with program() as prog:
            fr = declare(int)
            n = declare(int)
            I = declare(fixed)
            Q = declare(fixed)
            I_stream = declare_stream()
            Q_stream = declare_stream()

            with for_(n,0,n<avgs,n+1):
                with for_(fr,fmin,fr<fmax,fr+df):
                    update_frequency("resonator",59.8e6)
                    update_frequency("qubit",400e6-26.4e6)
                    wait(int(reset_time//4), "qubit")# wait for the qubit to relax, several T1s
                    if pi==True:
                        play("saturation","qubit",duration=40//4)
                    elif pi==False:
                        pass
                    update_frequency("qubit",fr)
                    play("saturation", "qubit",duration =1e3//4)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(0.6),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1")
                           )
                    save(I,I_stream)
                    save(Q,Q_stream)

            with stream_processing():
                I_stream.buffer(len(fvec)).average().save("I")
                Q_stream.buffer(len(fvec)).average().save("Q")
                
        return prog,fvec

    def ef_spec_V(self, *args):
        fmin = -473e6
        fmax = -468e6
        df = 0.05e6
        reset_time = 10e3
        fvec = np.arange(fmin, fmax, df)

        amin = -1.0
        amax = 1.0
        da = 0.05
        avec = np.arange(amin, amax + da / 2, da)
        reset_time = 15e3

        avgs = 5e3
        with program() as prog:
            fr = declare(int)
            n = declare(int)
            a=declare(fixed)
            I = declare(fixed)
            Q = declare(fixed)
            I_stream = declare_stream()
            Q_stream = declare_stream()

            with for_(n, 0, n < avgs, n + 1):
                with for_(fr, fmin, fr < fmax, fr + df):
                    with for_(a, amin, a < amax+da/2, a + da):

                        update_frequency("resonator", 59.8e6)
                        update_frequency("qubit", -24.5e6)
                        wait(int(reset_time // 4), "qubit")  # wait for the qubit to relax, several T1s
                      #  play("saturation", "qubit", duration=40 // 4)
                        update_frequency("qubit", fr)
                        play("saturation", "qubit", duration=1e3 // 4)
                        align("Vsource", "resonator","qubit")
                        play("CW" * amp(a), "Vsource", duration=int(5e3 // 4))  # 0.8 V for amp=1
                        reset_phase('resonator')
                        measure("readout" * amp(0.6), 'resonator', None,
                                demod.full('integW1', I, "out1"),
                                demod.full('integW2', Q, "out1")
                                )
                        save(I, I_stream)
                        save(Q, Q_stream)

            with stream_processing():
                I_stream.buffer(len(avec)).buffer(len(fvec)).average().save("I")
                Q_stream.buffer(len(avec)).buffer(len(fvec)).average().save("Q")

        return prog, fvec, avec

    def histogram(self,*args):
        T1=1e3
        reset_time = 5 * T1
        avgs = 100
        nvec = np.arange(0, avgs, 1)

        amin = 0.4
        amax = 0.6
        da = 0.05
        avec = np.arange(amin, amax + da / 2, da)

        with program() as prog:
            ##############################
            # declare real-time variables:
            ##############################

            n = declare(int)  #
            m = declare(int)  # Averaging
            a = declare(fixed)
            f = declare(int)

            I = declare(fixed)
            Q = declare(fixed)
            I1 = declare(fixed)
            Q1 = declare(fixed)
            I2 = declare(fixed)
            Q2 = declare(fixed)

            I_st = declare_stream()
            Q_st = declare_stream()

            #    update_frequency('resonator',10.2e6)
            update_frequency('qubit', 400e6-26.4e6)

            with for_(m, 0, m < 10000, m + 1):
                with for_(n, 0, n < avgs, n + 1):
                    with for_(a, amin, a < amax + da / 2, a + da):
                        update_frequency('resonator', 59.8e6)
                        #  """Just readout without playing anything"""
                        wait(int(reset_time // 4), "qubit")
                        align("qubit", "resonator")
                        reset_phase('resonator')
                        measure("readout" * amp(a), 'resonator', None,
                                demod.full('integW1', I1, "out1"),
                                demod.full('integW2', Q1, "out1"))

                        #             assign(I, I1-Q2)
                        #             assign(Q, I2+Q1)
                        save(I1, I_st)
                        save(Q1, Q_st)

                        #   """Play a ge pi pulse and then readout"""
                        wait(int(reset_time // 4), "qubit")
                        play("saturation", "qubit", duration=40//4)
                        align("qubit", "resonator")
                        reset_phase('resonator')
                        measure("readout" * amp(a), 'resonator', None,
                                demod.full('integW1', I1, "out1"),
                                demod.full('integW2', Q1, "out1"))

                        #             assign(I, I1-Q2)
                        #             assign(Q, I2+Q1)
                        save(I1, I_st)
                        save(Q1, Q_st)

            with stream_processing():
                I_st.buffer(2 * len(avec)).buffer(len(nvec)).average().save('I')
                Q_st.buffer(2 * len(avec)).buffer(len(nvec)).average().save('Q')
        return prog,avec, avgs

    def histogram_f(self,*args):
        T1=1e3
        reset_time = 5 * T1
        avgs = 1e4
        nvec = np.arange(0, avgs, 1)

        amin = 0.5
        amax = 0.7
        da = 0.05
        avec = np.arange(amin, amax + da / 2, da)

        fmin = 59.5e6
        fmax = 60.5e6
        df = 0.1e6
        fvec = np.arange(fmin, fmax, df)

        with program() as prog:
            ##############################
            # declare real-time variables:
            ##############################

            n = declare(int)  # Averaging
            a = declare(fixed)
            f = declare(int)

            I = declare(fixed)
            Q = declare(fixed)
            I1 = declare(fixed)
            Q1 = declare(fixed)
            I2 = declare(fixed)
            Q2 = declare(fixed)

            I_st = declare_stream()
            Q_st = declare_stream()

            #    update_frequency('resonator',10.2e6)
            update_frequency('qubit', 400e6-26.4e6)

            with for_(n, 0, n < avgs, n + 1):
                with for_(f, fmin, f < fmax , f + df):
                    update_frequency('resonator', f)
                    #  """Just readout without playing anything"""
                    wait(int(reset_time // 4), "qubit")
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout" * amp(0.6), 'resonator', None,
                            demod.full('integW1', I1, "out1"),
                            demod.full('integW2', Q1, "out1"))

                    #             assign(I, I1-Q2)
                    #             assign(Q, I2+Q1)
                    save(I1, I_st)
                    save(Q1, Q_st)

                    #   """Play a ge pi pulse and then readout"""
                    wait(int(reset_time // 4), "qubit")
                    play("saturation", "qubit", duration=40//4)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout" * amp(0.6), 'resonator', None,
                            demod.full('integW1', I1, "out1"),
                            demod.full('integW2', Q1, "out1"))

                    #             assign(I, I1-Q2)
                    #             assign(Q, I2+Q1)
                    save(I1, I_st)
                    save(Q1, Q_st)

            with stream_processing():
                I_st.buffer(2 * len(fvec)).buffer(len(nvec)).save('I')
                Q_st.buffer(2 * len(fvec)).buffer(len(nvec)).save('Q')
        return prog,fvec

    def histogram_fid(self,Is, Qs, numbins=100):

        a_num = len(Is[0])

        colors = ['r', 'b', 'g', 'orange']
        labels = ['g', 'e', 'f', 'h']
        titles = ['-I', '-Q']

        fig = plt.figure(figsize=(12, 16))

        ax = fig.add_subplot(421)
        x0g, y0g = np.mean(Is[0]), np.mean(Qs[0])
        x0e, y0e = np.mean(Is[1]), np.mean(Qs[1])
        phi = np.arctan((y0e - y0g) / (x0e - x0g))
        for ii in range(2):
            ax.plot(Is[ii], Qs[ii], '.', color=colors[ii], alpha=0.85)
        ax.plot(x0g, y0g, 'v', color='black')
        ax.plot(x0e, y0e, '^', color='black')

        ax.set_xlabel('I (V)')
        ax.set_ylabel('Q (V)')

        ax = fig.add_subplot(422)

        ax.plot(x0g, y0g, 'v', color='black')
        ax.plot(x0e, y0e, '^', color='black')
        ax.set_xlabel('I')
        ax.set_ylabel('Q')

        Isrot = [Is[ii] * np.cos(phi) + Qs[ii] * np.sin(phi) for ii in range(2)]
        Qsrot = [-Is[ii] * np.sin(phi) + Qs[ii] * np.cos(phi) for ii in range(2)]

        ax = fig.add_subplot(423, title='rotated')
        Is, Qs = Isrot, Qsrot

        x0g, y0g = np.mean(Is[0]), np.mean(Qs[0])
        x0e, y0e = np.mean(Is[1]), np.mean(Qs[1])
        phi = np.arctan((y0e - y0g) / (x0e - x0g))
        for ii in range(2):
            ax.plot(Is[ii], Qs[ii], '.', color=colors[ii], alpha=0.85)
        ax.plot(x0g, y0g, 'v', color='black')
        ax.plot(x0e, y0e, '^', color='black')

        ax.set_xlabel('I (V)')
        ax.set_ylabel('Q (V)')

        ax = fig.add_subplot(424)

        ax.plot(x0g, y0g, 'v', color='black')
        ax.plot(x0e, y0e, '^', color='black')
        ax.set_xlabel('I')
        ax.set_ylabel('Q')

        ax = fig.add_subplot(4, 2, 5, title='I')
        ax.hist(Is[0], bins=numbins, alpha=0.75, color=colors[0])
        ax.hist(Is[1], bins=numbins, alpha=0.75, color=colors[1])
        ax.set_xlabel('I' + '(V)')
        ax.set_ylabel('Number')
        ax.legend()

        ax = fig.add_subplot(4, 2, 6, title='Q')
        ax.hist(Qs[0], bins=numbins, alpha=0.75, color=colors[0])
        ax.hist(Qs[1], bins=numbins, alpha=0.75, color=colors[1])

        ax.set_xlabel('Q' + '(V)')
        ax.set_ylabel('Number')
        ax.legend()

        sshg, ssbinsg = np.histogram(Is[0], bins=numbins, range=[min(Is[0]), max(Is[0])])
        sshe, ssbinse = np.histogram(Is[1], bins=numbins, range=[min(Is[0]), max(Is[0])])
        fid = np.abs(((np.cumsum(sshg) - np.cumsum(sshe)) / sshg.sum())).max()

        returnfid = fid
        print("Single shot readout fidility from channel I", " = ", fid)
        print('---------------------------')

        ax = fig.add_subplot(4, 2, 7)
        ax.plot(ssbinse[:-1], np.cumsum(sshg) / sshg.sum(), color=colors[0])
        ax.plot(ssbinse[:-1], np.cumsum(sshe) / sshg.sum(), color=colors[1])
        ax.plot(ssbinse[:-1], np.abs(np.cumsum(sshe) - np.cumsum(sshg)) / sshg.sum(), color='k')

        sshg, ssbinsg = np.histogram(Qs[0], bins=numbins, range=[min(Qs[0]), max(Qs[0])])
        sshe, ssbinse = np.histogram(Qs[1], bins=numbins, range=[min(Qs[0]), max(Qs[0])])
        fid = np.abs(((np.cumsum(sshg) - np.cumsum(sshe)) / sshg.sum())).max()

        print("Single shot readout fidility from channel Q", i, " = ", fid)
        print('---------------------------')

        ax = fig.add_subplot(4, 2, 8)
        ax.plot(ssbinse[:-1], np.cumsum(sshg) / sshg.sum(), color=colors[0])
        ax.plot(ssbinse[:-1], np.cumsum(sshe) / sshg.sum(), color=colors[1])
        ax.plot(ssbinse[:-1], np.abs(np.cumsum(sshe) - np.cumsum(sshg)) / sshg.sum(), color='k')

        fig.tight_layout()
        plt.show()

        return (returnfid)

    def histogram_simp_fid(self, Is, Qs, ax, numbins=100):

        a_num = len(Is[0])

        colors = ['r', 'b', 'g', 'orange']
        labels = ['g', 'e', 'f', 'h']
        titles = ['-I', '-Q']

        x0g, y0g = np.mean(Is[0]), np.mean(Qs[0])
        x0e, y0e = np.mean(Is[1]), np.mean(Qs[1])
        phi = np.arctan((y0e - y0g) / (x0e - x0g))

        Isrot = [Is[ii] * np.cos(phi) + Qs[ii] * np.sin(phi) for ii in range(2)]
        Qsrot = [-Is[ii] * np.sin(phi) + Qs[ii] * np.cos(phi) for ii in range(2)]

        Is, Qs = Isrot, Qsrot

        x0g, y0g = np.mean(Is[0]), np.mean(Qs[0])
        x0e, y0e = np.mean(Is[1]), np.mean(Qs[1])
        phi = np.arctan((y0e - y0g) / (x0e - x0g))

        sshg, ssbinsg = np.histogram(Is[0], bins=numbins, range=[min(Is[0]), max(Is[0])])
        sshe, ssbinse = np.histogram(Is[1], bins=numbins, range=[min(Is[0]), max(Is[0])])
        fid = np.abs(((np.cumsum(sshg) - np.cumsum(sshe)) / sshg.sum())).max()

        ax.plot(ssbinse[:-1], np.cumsum(sshg) / sshg.sum(), color=colors[0])
        ax.plot(ssbinse[:-1], np.cumsum(sshe) / sshg.sum(), color=colors[1])
        ax.plot(ssbinse[:-1], np.abs(np.cumsum(sshe) - np.cumsum(sshg)) / sshg.sum(), color='k')

        return (fid)

    def ef_rabi2d(self,*args):
        dt = 2
        T_min = 1
        T_max = 300//4
        times = np.arange(T_min, T_max, dt)*4

        T1=1e3

        fmin = 200e6
        fmax = 250e6
        df = 0.1e6
        fvec=np.arange(fmin,fmax,df)
        reset_time = 5*T1

        avgs = 50000
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
            update_frequency("resonator",59.8e6)
            with for_(n, 0, n < avgs, n + 1):
                with for_(f, fmin,f<fmax, f+df):
                    with for_(t, T_min, t < T_max, t + dt):
                        wait(int(reset_time//4), "qubit")
                        update_frequency("qubit", 400e6-24.5e6)
                        play("saturation","qubit",duration=40//4)
                        update_frequency("qubit", f)
                        play("saturation", "qubit",duration = t)
                        align("qubit", "resonator")
                        reset_phase('resonator')
                        measure("readout"*amp(0.6),'resonator',None,
                                demod.full('integW1',I,"out1"),
                                demod.full('integW2',Q,"out1")
                               )
                        save(I,I_st)
                        save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)).buffer(len(fvec)).average().save('I')
                Q_st.buffer(len(times)).buffer(len(fvec)).average().save('Q')

        return prog,fvec +3.7e9,times

    def ef_rabi1d(self,*args):
        dt = 1
        T_min = 1
        T_max = 200//4
        times = np.arange(T_min, T_max, dt)*4

        T1=10e3
        reset_time = 5*T1

        avgs = 1e6
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
            update_frequency("resonator",59.8e6)
            with for_(n, 0, n < avgs, n + 1):
                with for_(t, T_min, t < T_max, t + dt):
                    wait(int(reset_time//4), "qubit")
                    update_frequency("qubit", -24.5e6)
                    play("saturation","qubit",duration=40//4)
                    update_frequency("qubit", -470e6)
                    play("saturation"*amp(0.7), "qubit",duration=t)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(0.6),'resonator',None,
                            demod.full('integW1',I,"out1"),
                            demod.full('integW2',Q,"out1")
                           )
                    save(I,I_st)
                    save(Q,Q_st)

                """Play a ge pi pulse and then readout"""
                wait(int(reset_time // 4), "qubit")
                update_frequency("qubit", -24.5e6)
                play("saturation", "qubit", duration=40 // 4)
                update_frequency("qubit", -529.8e6)
                play("saturation"*amp(0.7), "qubit", duration=380 // 4)
                align('qubit', 'resonator')
                reset_phase('resonator')
                measure("readout" * amp(0.6), 'resonator', None,
                        demod.full('integW1', I, "out1"),
                        demod.full('integW2', Q, "out1"))
                save(I, I_st)
                save(Q, Q_st)

                """Just readout without playing anything"""
                wait(int(reset_time // 4), "qubit")
                update_frequency("qubit", -24.5e6)
                play("saturation", "qubit", duration=40 // 4)
                align('qubit', 'resonator')
                reset_phase('resonator')
                measure("readout" * amp(0.6), 'resonator', None,
                        demod.full('integW1', I, "out1"),
                        demod.full('integW2', Q, "out1"))
                save(I, I_st)
                save(Q, Q_st)

            with stream_processing():
                I_st.buffer(len(times)+2).average().save('I')
                Q_st.buffer(len(times)+2).average().save('Q')

        return prog,times

    def ef_T1(self,*args):

        dt = 20
        T_min = 1
        T_max = 1000
        times = np.arange(T_min, T_max, dt)*4
        T1=10.5e3
        reset_time = 5*T1

        avgs = 1000
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
            update_frequency("resonator",59.8e6)
            with for_(n, 0, n < avgs, n + 1):
                with for_(t, T_min, t < T_max, t + dt):
                    wait(int(reset_time//4), "qubit")
                    update_frequency("qubit",-24.5e6)
                    play("saturation", "qubit",duration=40//4) # pi pulse with saturation
                    update_frequency("qubit",-529.8e6)
                    play("saturation","qubit",duration=200//4)
                    wait(t)
                    align("qubit", "resonator")
                    reset_phase('resonator')
                    measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))
                    save(I,I_st)
                    save(Q,Q_st)

                wait(int(reset_time//4), "qubit")
                play("saturation", "qubit", duration=40// 4)  # pi pulse with saturation
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None, #amp=0.45, f=10.2 looks best for low power (resonator punched out) readout
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1")
                               )
                save(I,I_st)
                save(Q,Q_st)

                wait(int(reset_time//4), "qubit")
                align("qubit", "resonator")
                reset_phase('resonator')
                measure("readout"*amp(0.6),'resonator',None,
                        demod.full('integW1',I,"out1"),
                        demod.full('integW2',Q,"out1"))

                save(I,I_st)
                save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)+2).average().save('I')
                Q_st.buffer(len(times)+2).average().save('Q')

        return prog,times

    def T1rho(self,*args):
        dt = 50//4
        T_min = 1
        T_max = 1000//4
        times = np.arange(T_min, T_max, dt)*4

        T1=1e3

        amin = 0
        amax = 1.0
        da = 0.01
        avec=np.arange(amin,amax,da)
        reset_time = 5*T1

        avgs = 50000
        with program() as prog:

            ##############################
            # declare real-time variables:
            ##############################

            n = declare(int)      # Averaging
            i = declare(int)      # Amplitudes
            t = declare(int)     #array of time delays
            a = declare(fixed)
            I = declare(fixed)
            Q = declare(fixed)

            I_st = declare_stream()
            Q_st = declare_stream()


            ###############
            # the sequence:
            ###############
            update_frequency("resonator",59.8e6)
            update_frequency("qubit",400e6-26.4e6)

            with for_(n, 0, n < avgs, n + 1):
                with for_(a, amin,a<amax+da/2, a+da):
                    with for_(t, T_min, t < T_max, t + dt):
                        wait(int(reset_time//4), "qubit")
                        play("saturation","qubit",duration=40//8)
                        play("saturation"*amp(a), "qubit",duration = t)
                        align("qubit", "resonator")
                        play("saturation","qubit",duration=40//8)
                        reset_phase('resonator')
                        measure("readout"*amp(0.6),'resonator',None,
                                demod.full('integW1',I,"out1"),
                                demod.full('integW2',Q,"out1")
                               )
                        save(I,I_st)
                        save(Q,Q_st)

            with stream_processing():
                I_st.buffer(len(times)).buffer(len(avec)).average().save('I')
                Q_st.buffer(len(times)).buffer(len(avec)).average().save('Q')

        return prog,avec,times
