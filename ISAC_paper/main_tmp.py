import numpy as np
import matlab.engine
import matlab


BEAM_PARAM = {

    "NUM_USER" : 4,
    "NUM_TARGET" : 2,
    "NUM_ANTENNA" : 12,
    "USER" : np.array([[-300, -300],[-70, -400], [70, -400], [300, -300]]),
    "TARGET" : np.array([[-100, 100], [100, 100]]),
    "UAV" : np.array([0, 0]),
    "SENSING_THRESHOLD_db" : -22,
    "UAV_Z" : 30,
    "V_MAX" : 30,
    "P_MAX" : 0.5,
    "CHANNEL_GAIN" : 10 ** (-6),
    "NOISE_POWER" : 10 ** (-14),
    "SCALING" : 10 ** (3)
}

BEAM_PARAM["SENSING_THRESHOLD"] = 10 ** (0.1 * BEAM_PARAM["SENSING_THRESHOLD_db"]) * 10 ** (-3)
BEAM_PARAM["NOISE_POWER_SCALING"] = BEAM_PARAM["NOISE_POWER"] * BEAM_PARAM["SCALING"] ** (2)
BEAM_PARAM["SENSING_THRESHOLD_SCALING"] = BEAM_PARAM["SENSING_THRESHOLD"] * BEAM_PARAM["SCALING"] ** (2)



if __name__ == "__main__":

    eng = matlab.engine.start_matlab()
    eng.addpath(r'C:\\Users\\Wireless\\Desktop\\개인연구\\other\\other_other\\seunghyeon_shin\\matlab_part')

    matlab_param = {key: matlab.double(value.tolist()) if isinstance(value, np.ndarray) else matlab.double([value]) for key, value in BEAM_PARAM.items()}
    matlab_param_struct = eng.struct(matlab_param)

    W_opt, R_opt = eng.ISAC_paper_BEAMFORMING_python(matlab_param_struct, nargout=2)
    W_opt = np.array(W_opt)
    R_opt = np.array(R_opt)

    eng.quit()

        
