import numpy as np
import json
import glob
import matplotlib.pyplot as plt
from scipy import stats
import math
from scipy.signal import butter, filtfilt
from scipy.io import savemat
import argparse
import os

def read_data(data_id):
    with open(f"./wam_data/state_{data_id}.json") as f:
        state_data = json.load(f)

    with open(f"./wam_data/gravity_{data_id}.json") as f:
        gravity_data = json.load(f)

    return state_data, gravity_data

def joint_traj_eval(t, f=0.1, A=0.6):
    f *= 2 * math.pi
    q = A * np.array([0, np.sin(f*t), 0, np.sin(f*t)], dtype=object)
    qdot = A * np.array([0, f * np.cos(f*t), 0, f * np.cos(f*t)], dtype=object)
    qdotdot = A * np.array([0, -f**2 * np.sin(f*t), 0, -f**2 * np.sin(f*t)], dtype=object)
    return q, qdot, qdotdot

def joint_traj_excite(t, T=10):
    t = np.asarray(t)  # ensure t is ndarray
    L = 5
    x = [
        0.707990384021940, -0.0150146020506483, 0.000280191811242826, -0.150424710183653, -0.542831263598882,
        0.0285165790096300, -0.0379420078788927, -0.000370896589483716, -0.126234977048364, 0.110684006942012,
        0.0679857584856881, -0.000160746280917686, 3.59730863728628e-05, -0.608666198402881, 0.540805213111737,
        -0.0433137281308812, -0.000693832290851609, -0.000400179294006102, 0.466199557165230, -0.363779259613263
    ]
    A2 = np.array(x[0:L])[:, np.newaxis]
    B2 = np.array(x[L:2*L])[:, np.newaxis]
    A4 = np.array(x[2*L:3*L])[:, np.newaxis]
    B4 = np.array(x[3*L:4*L])[:, np.newaxis]

    omega = 2 * np.pi / T
    k = np.arange(1, L+1)[:, np.newaxis]  # shape (L, 1)

    omega_kt = omega * k * t  # shape (L, N)

    sin_term = np.sin(omega_kt)
    cos_term = np.cos(omega_kt)

    # Positions
    theta2 = np.sum((A2 * sin_term - B2 * cos_term) / (omega * k), axis=0)
    theta4 = np.sum((A4 * sin_term - B4 * cos_term) / (omega * k), axis=0)

    # Velocities
    thetad2 = np.sum((A2 * cos_term + B2 * sin_term), axis=0)
    thetad4 = np.sum((A4 * cos_term + B4 * sin_term), axis=0)

    # Accelerations
    thetadd2 = np.sum((-A2 * sin_term + B2 * cos_term) * (omega * k), axis=0)
    thetadd4 = np.sum((-A4 * sin_term + B4 * cos_term) * (omega * k), axis=0)

    factor = 1.0
    q = np.vstack([np.zeros_like(t), theta2, np.zeros_like(t), theta4]) * factor
    qdot = np.vstack([np.zeros_like(t), thetad2, np.zeros_like(t), thetad4]) * factor
    qdotdot = np.vstack([np.zeros_like(t), thetadd2, np.zeros_like(t), thetadd4]) * factor

    return q, qdot, qdotdot

def butter_lowpass_filter(data, cutoff=2, nyq=0.5*300, order=2):
    normal_cutoff = cutoff/nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def FFD_lowpass_filter(qo, fs=300, kt=3):
    """
    Apply Frequency-domain Filter and Differential (FFD) to obtain time derivative of joint angles without high-frequency noise.
    
    Parameters:
    qo : numpy array
        Average joint angles.
    fs : float
        Sampling frequency.
    kt : float
        Threshold frequency for the filter.
    
    Returns:
    q_dot_FFD : numpy array
        Joint velocity.
    q_ddot_FFD : numpy array
        Joint acceleration.
    """
    # Perform Discrete Fourier Transform (DFT)
    F = np.fft.fft(qo)
    freqs = np.fft.fftfreq(len(qo), 1/fs)
    # print(freqs)

    # Apply frequency-domain filter
    Ff = F.copy()
    Ff[np.abs(freqs) > kt] = 0

    # Perform Inverse DFT to get noise-free data
    qof = np.fft.ifft(Ff).real

    # Compute cyclic frequencies
    omega = 2 * np.pi * freqs

    # Calculate first time-derivative (velocity)
    q_dot_FFD = np.fft.ifft(Ff * (1j * omega)).real

    # Calculate second time-derivative (acceleration)
    q_ddot_FFD = np.fft.ifft(Ff * (1j * omega) ** 2).real

    return qof, q_dot_FFD, q_ddot_FFD

def central_difference(t, y):
    if len(t) != len(y):
        raise ValueError("x and y must have the same length")

    dydt = np.zeros_like(y)

    # Central difference for the interior points
    for i in range(1, len(t)-1):
        dydt[i] = (y[i+1] - y[i-1]) / (t[i+1] - t[i-1])

    return dydt

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("data_id", type=int)
    args = parser.parse_args()

    data_id = args.data_id
    state_data, gravity_data = read_data(data_id)

    time, pos2, vel2, torq2, grav2 = [], [], [], [], []
    pos4, vel4, torq4, grav4 = [], [], [], []

    for s, g in zip(state_data, gravity_data):
        time.append(s['time'])

        pos2.append(s['position'][1])
        vel2.append(s['velocity'][1])
        torq2.append(s['effort'][1])
        grav2.append(g['gravity'][1])

        pos4.append(s['position'][3])
        vel4.append(s['velocity'][3])
        torq4.append(s['effort'][3])
        grav4.append(g['gravity'][3])

    time = np.array(time)
    time -= time[0]

    pos = np.array([pos2, pos4])
    vel = np.array([vel2, vel4])
    torq = np.array([torq2, torq4])
    grav = np.array([grav2, grav4])

    pos2_freq = np.fft.fft(pos[0])
    pos4_freq = np.fft.fft(pos[1])

    frequency = len(time)/time[-1]
    print("frequency: ", frequency)

    diff_time = np.diff(time)
##    print("time diff: ", diff_time)


    vel2_prcd, acc2_prcd, _ = FFD_lowpass_filter(vel[0])
    vel4_prcd, acc4_prcd, _ = FFD_lowpass_filter(vel[1])

    torq2_prcd, _, _ = FFD_lowpass_filter(torq[0])
    torq4_prcd, _, _ = FFD_lowpass_filter(torq[1])

    grav2_prcd, _, _ = FFD_lowpass_filter(grav[0])
    grav4_prcd, _, _ = FFD_lowpass_filter(grav[1])

    pos_des, vel_des, acc_des = joint_traj_excite(time)
    
    plt.subplot(1, 2, 1)
    plt.plot(time, pos2_freq.real)
    plt.title('joint 2')
    plt.subplot(1, 2, 2)
    plt.plot(time, pos4_freq.real)
    plt.title('joint 4')
    plt.suptitle('position')
    plt.show()
    
    plt.subplot(1, 2, 1)
    plt.plot(time, pos[0])
    plt.plot(time, pos_des[1])
##    plt.plot(time, pos2_prcd)
    plt.legend(['collected', 'desired'])
    plt.title('joint 2')
    plt.subplot(1, 2, 2)
    plt.plot(time, pos[1])
    plt.plot(time, pos_des[3])
##    plt.plot(time, pos4_prcd)
    plt.legend(['collected', 'desired'])
    plt.title('joint 4')
    plt.suptitle('position')
    plt.show()

    plt.subplot(1, 2, 1)
    plt.plot(time, vel[0])
    plt.plot(time, vel_des[1])
    plt.plot(time, vel2_prcd)
    plt.legend(['collected', 'desired', 'filtered'])
    plt.title('joint 2')
    plt.subplot(1, 2, 2)
    plt.plot(time, vel[1])
    plt.plot(time, vel_des[3])
    plt.plot(time, vel4_prcd)
    plt.legend(['collected', 'desired', 'filtered'])
    plt.title('joint 4')
    plt.suptitle('velocity')
    plt.show()

    plt.subplot(1, 2, 1)
##    plt.plot(time, acc2)
    plt.plot(time, acc_des[1])
    plt.plot(time, acc2_prcd)
    plt.legend(['desired', 'filtered'])
    plt.title('joint 2')
    plt.subplot(1, 2, 2)
##    plt.plot(time, acc4)
    plt.plot(time, acc_des[3])
    plt.plot(time, acc4_prcd)
    plt.legend(['desired', 'filtered'])
    plt.title('joint 4')
    plt.suptitle('acceleration')
    plt.show()

    plt.subplot(1, 2, 1)
    plt.plot(time, torq[0])
    plt.plot(time, torq2_prcd)
    plt.legend(['collected', 'filtered'])
    plt.title('joint 2')
    plt.subplot(1, 2, 2)
    plt.plot(time, torq[1])
    plt.plot(time, torq4_prcd)
    plt.legend(['collected', 'filtered'])
    plt.title('joint 4')
    plt.suptitle('torque')
    plt.show()

    plt.subplot(1, 2, 1)
    plt.plot(time, grav[0])
    plt.plot(time, grav2_prcd)
    plt.legend(['collected', 'filtered'])
    plt.title('joint 2')
    plt.subplot(1, 2, 2)
    plt.plot(time, grav[1])
    plt.plot(time, grav4_prcd)
    plt.legend(['collected', 'filtered'])
    plt.title('joint 4')
    plt.suptitle('gravity')
    plt.show()

    processed_data = [pos[0], pos[1], vel2_prcd, vel4_prcd, acc2_prcd, acc4_prcd, torq2_prcd, torq4_prcd, grav2_prcd, grav4_prcd]
    processed_data = np.array(processed_data)
    output_dir = "./processed_data"
    os.makedirs(output_dir, exist_ok=True)
    savemat(f"{output_dir}/processed_data_{data_id}.mat", {'processed_data': processed_data})


    




    


