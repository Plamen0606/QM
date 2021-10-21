import numpy as np

######################
# AUXILIARY FUNCTIONS:
######################


def gauss(amplitude, mu, sigma, length):
    t = np.linspace(-length / 2, length / 2, length)
    gauss_wave = amplitude * np.exp(-((t - mu) ** 2) / (2 * sigma ** 2))
    return [float(x) for x in gauss_wave]


def IQ_imbalance(g, phi):
    c = np.cos(phi)
    s = np.sin(phi)
    N = 1 / ((1-g**2)*(2*c**2-1))
    return [float(N * x) for x in [(1-g)*c, (1+g)*s, (1-g)*s, (1+g)*c]]


################
# CONFIGURATION:
################

long_readout_len = 10000
readout_len = 5000

qubit_IF = 150e6
rr_IF =  30e6  # rr -> readout resonator

qubit_LO = 6.345e9

gauss_len=200

rr_freq = 7.26e9
rr_LO = rr_freq - rr_IF

config = {

    'version': 1,

    'controllers': {
        'con1': {
            'type': 'opx1',
            'analog_outputs': {
				1: {'offset': 0.0},  #Voltage source I
				2: {'offset': 0.0},  #Voltage source Q
                5: {'offset': -0.025},  # qubit I
                6: {'offset': -0.008},  # qubit Q
                3: {'offset': 0.0},  # RR I
                4: {'offset': 0.0},  # RR Q
            },
            'digital_outputs': {
				2: {},
			},
            'analog_inputs': {
                1: {'offset': 0.17885652719067774,"gain_db":2},  # I
                2: {'offset': 0.0}   # Q
            }
        }
    },

    'elements': {

        'qubit': {
            'mixInputs': {
                'I': ('con1', 5),
                'Q': ('con1', 6),
                'lo_frequency': qubit_LO,
                'mixer': 'mixer_qubit'
            },
            'intermediate_frequency': qubit_IF,
            'operations': {
                'CW': 'CW',
                'saturation': 'saturation_pulse',
                'gaussian': 'gaussian_pulse',
                'pi': 'pi_pulse',
				'marker': 'marker_pulse',
            },
			'digitalInputs': {
                "Switch": {
                    "port": ("con1", 2),
                    "delay": 0,
                    "buffer": 0,
                },
            },

        },

        'resonator': {
            'mixInputs': {
                'I': ('con1', 3),
                'Q': ('con1', 4),
                'lo_frequency': rr_LO,
                'mixer': 'mixer_RR'
            },
            'intermediate_frequency': rr_IF,
            'operations': {
                'CW': 'CW',
                'long_readout': 'long_readout_pulse',
                'readout': 'readout_pulse',
            },
            "outputs": {
                'out1': ('con1', 1),
                'out2': ('con1', 2)
            },
            'time_of_flight': 200,
            'smearing': 0
        },
		'Vsource': {
            'mixInputs': {
                'I': ('con1', 1),
                'Q': ('con1', 2),
                'lo_frequency': 0,
            },
            'intermediate_frequency': 0,
            'operations': {
                'CW': 'CW',
                'long_readout': 'long_readout_pulse',
                'readout': 'readout_pulse',
            },
            "outputs": {
                'out1': ('con1', 1),
                'out2': ('con1', 2)
            },
            'time_of_flight': 200,
            'smearing': 0
        },
		
    },

    "pulses": {

        "CW": {
            'operation': 'control',
            'length': 60000,
            'waveforms': {
                'I': 'const_wf',
                'Q': 'zero_wf'
            }
        },

		"marker_pulse": {
            'operation': 'control',
            'length': 20,
            'waveforms': {
                'I': 'zero_wf',
                'Q': 'zero_wf'			
            },
			'digital_marker': 'ON',

        },
		
        "saturation_pulse": {
            'operation': 'control',
            'length': 50000,  # several T1s
            'waveforms': {
                'I': 'saturation_wf',
                'Q': 'zero_wf'
            }
        },

        "gaussian_pulse": {
            'operation': 'control',
            'length': gauss_len,
            'waveforms': {
                'I': 'gauss_wf',
                'Q': 'zero_wf'
            }
        },

        'pi_pulse': {
            'operation': 'control',
            'length': 6000,
            'waveforms': {
                'I': 'pi_wf',
                'Q': 'zero_wf'
            }
        },

        'long_readout_pulse': {
            'operation': 'measurement',
            'length': long_readout_len,
            'waveforms': {
                'I': 'long_readout_wf',
                'Q': 'zero_wf'
            },
            'integration_weights': {
                'long_integW_cos': 'long_integW_cos',
                'long_integW_sin': 'long_integW_sin',
            },
			"digital_marker":"ON",
        },

        'readout_pulse': {
            'operation': 'measurement',
            'length': readout_len,
            'waveforms': {
                'I': 'readout_wf',
                'Q': 'zero_wf'
            },
            'integration_weights': {
                'integW1': 'integW1',
                'integW2': 'integW2',
                'optW1': 'optW1',
                'optW2': 'optW2'
            },
			"digital_marker":"ON",
        },
		'long_empty_readout_pulse':{
            'operation': 'measurement',
            'length': long_readout_len,
            'waveforms': {
                'I': 'zero_wf',
                'Q': 'zero_wf'
            },
            'integration_weights': {
                'long_integW_cos': 'long_integW_cos',
                'long_integW_sin': 'long_integW_sin',
            },
			"digital_marker":"ON",
        },

    },

    'waveforms': {

        'const_wf': {
            'type': 'constant',
            'sample': 0.4
        },

        'zero_wf': {
            'type': 'constant',
            'sample': 0.0
        },

        'saturation_wf': {
            'type': 'constant',
            'sample': 0.49
        },

        'gauss_wf': {
            'type': 'arbitrary',
            'samples': gauss(0.49, 0.0, gauss_len//4, gauss_len)
        },

        'pi_wf': {
            'type': 'arbitrary',
            'samples': gauss(0.3, 0.0, 6.0, 6000)
        },

        'long_readout_wf': {
            'type': 'constant',
            'sample': 0.32
        },

        'readout_wf': {
            'type': 'constant',
            'sample': 0.34
        },
    },

    'digital_waveforms': {
        'ON': {
            'samples': [(1, 0)]
        }
    },

    'integration_weights': {

        'long_integW_cos': {
            'cosine': [1.0] * int(long_readout_len / 4),
            'sine': [0.0] * int(long_readout_len / 4)
        },

        'long_integW_sin': {
            'cosine': [0.0] * int(long_readout_len / 4),
            'sine': [1.0] * int(long_readout_len / 4)
        },

        'integW1': {
            'cosine': [1.0] * int(readout_len / 4),
            'sine': [0.0] * int(readout_len / 4),
        },

        'integW2': {
            'cosine': [0.0] * int(readout_len / 4),
            'sine': [1.0] * int(readout_len / 4),
        },

        'optW1': {
            'cosine': [1.0] * int(readout_len / 4),
            'sine': [0.0] * int(readout_len / 4)
        },

        'optW2': {
            'cosine': [0.0] * int(readout_len / 4),
            'sine': [1.0] * int(readout_len / 4)
        },
    },

    'mixers': {
        'mixer_qubit': [
            {'intermediate_frequency': qubit_IF, 'lo_frequency': qubit_LO,
             'correction':  IQ_imbalance(-0.01, -0.02)},
        ],
        'mixer_RR': [
            {'intermediate_frequency': rr_IF, 'lo_frequency': rr_LO,
             'correction': IQ_imbalance(0.0, 0.0)}
        ],
    }
}
