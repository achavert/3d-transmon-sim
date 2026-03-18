import numpy as np
from scipy.integrate import RK45

class Transmon3D:
    #Constants
    kB = 1.380649e-23 # J K-1
    hbar = 6.626070151e-34 # J Hz-1

    X01 = np.array([[0, 1, 0],
                    [1, 0, 0],
                    [0, 0, 0]], dtype = complex)
    
    Y01 = np.array([[0, -1j, 0],
                    [1j, 0, 0],
                    [0, 0, 0]], dtype = complex)
    
    X12 = np.array([[0, 0, 0],
                    [0, 0, 1],
                    [0, 1, 0]], dtype = complex)
    
    Y12 = np.array([[0, 0, 0],
                    [0, 0, -1j],
                    [0, 1j, 0]], dtype = complex)

    def __init__(self,
                 freq : float = 4.3910e+09, #Hz
                 anh : float = -1.8100e+08, #Hz
                 lambda1 : float = 1.0,
                 lambda2 : float = 1.0,
                 dt : float = 0.5e-9, #s
                 **noise):
        
        self.freq = freq
        self.anh = anh
        self.lambda1 = lambda1
        self.lambda2 = lambda2
        self.noise = noise

        self.init_state = self.build_initial_state()
        self.current_state = self.init_state

        self.schedules = []
        self.time = 0
        self.total_time = 0
    
    def build_initial_state(self):
        if len(self.noise.values()) == 0:
            return np.array([[1, 0, 0],
                             [0, 0, 0],
                             [0, 0, 0]], 
                            dtype = complex)
        else :
            return NotImplementedError('Noise simulation not implemented yet')

    def add_new_schedule(self, sched : Schedule3D):
        self.schedules.append(sched)
        self.total_time += sched.dur
    
    def solve_schedule(self, sched):
        func = self.get_func(sched)
        solver = RK45(fun = func, 
                      t0 = self.time * self.dt, 
                      y0 = self.y0, 
                      t_bound = (self.time + sched.dur) * self.dt,
                      rtol = 1e-4,
                      atol = 1e-7)
        
        y_points = []
        while solver.status == 'running':
            self.time += 1
            solver.step() 
            y_points.append(solver.y.copy())
        y_points = np.array(y_points)
        return y_points



class Schedule3D:
    def __init__(self, drive_freq, amp, sch):
        self.drive_freq = drive_freq
        self.amp = amp
        self.sch = sch
        if len(sch[0]) == len(sch[1]):
            self.dur = len(sch[0])
        else:
            raise Exception('Schedule components do not match in length') 

    def get_components(self, step):
        return (self.sch(step), self.sch(step))
    
