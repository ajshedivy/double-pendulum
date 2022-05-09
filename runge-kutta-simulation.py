from pathlib import Path
from typing import List, Union
from xmlrpc.client import DateTime
from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from collections import deque
import argparse
import pandas as pd

G = 9.8  # acceleration due to gravity, in m/s^2
L1 = 1.0  # length of pendulum in m
L = L1 + L1  # maximal length of the combined pendulum
history_len = 500  # how many trajectory points to display

fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, L))
ax.set_aspect('equal')
ax.grid()

Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800)

line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], '.-', lw=1, ms=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)

def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    if i == 0:
        history_x.clear()
        history_y.clear()

    history_x.appendleft(thisx[2])
    history_y.appendleft(thisy[2])

    line.set_data(thisx, thisy)
    trace.set_data(history_x, history_y)
    time_text.set_text(time_template % (i*dt))
    return line, trace, time_text

def read_data(filename: Union[str, Path]) -> List[Union[str, float]]:
    data = list()
    with open(filename, 'r') as f:
        data = f.readlines()
        data = [line.rstrip().split() for line in data]

    return data

def get_aminations(data: List[str], gif: Union[str, Path], video: Union[str, Path]) -> None:
    global x1, x2, y1, y2, dt
    df = pd.DataFrame(data)
    _x1 = 0
    _y1 = 1
    _x2 = 2
    _y2 = 3
    _t = 4
    _theta1 = 5
    _theta2 = 6

    x1 = df[_x1].to_numpy().astype(float)
    y1 = df[_y1].to_numpy().astype(float)
    x2 = df[_x2].to_numpy().astype(float)
    y2 = df[_y2].to_numpy().astype(float)
    t = df[_t].to_numpy().astype(float)
    dt = 20 / len(x1)
    theta1 = df[_theta1].to_numpy().astype(float)
    theta2 = df[_theta2].to_numpy().astype(float)
    history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)

    ani = animation.FuncAnimation(
        fig, animate, len(x1), interval=dt*1000, blit=True)

    ani.save(gif)
    ani.save(video, writer=writer)

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', default='output.txt')
    parser.add_argument('--gif', default='simulation.gif')
    parser.add_argument('--video', default='simulation.mp4')
    return parser

def parse_args() -> argparse.Namespace:
    parser = get_parser()
    opts   = parser.parse_args()
    return opts


def main():
    opts          = parse_args()
    data_file     = opts.data
    gif_file      = opts.gif
    mp4_file      = opts.video

    data = read_data(data_file)

    get_aminations(data, gif_file, mp4_file)

if __name__ == '__main__':
    main()

    

