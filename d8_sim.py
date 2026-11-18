"""ME5250 Project 1, D8: Peaucellier-Lipkin linkage in MuJoCo (independent check).

Script version of d8_sim.ipynb. This file does not call the project FK code.

Tree (7 kept joints, one per moving body):
  world --O--> OA --A--> DA
  world --O--> OC --C--> CD
  world --E--> EB --B--> AB
                    \\--B--> BC
Cut joints (one per loop), closed with <connect>:
  A: end of AB == end of OA
  C: end of BC == end of OC
  D: end of DA == end of CD

Outputs: d8_mujoco_log.csv (settled measurements), d8_video.mp4 (top view).
Run:  .venv/bin/python d8_sim.py
"""
import numpy as np
import mujoco
import imageio


r, L, s = 60.0, 150.0, 60.0
MM = 1e-3



xml_one = f"""
<mujoco>
  <worldbody>
    <body name="OA" pos="0 0 0">
      <joint name="OA" type="hinge" axis="0 0 1"/>
      <geom type="capsule" fromto="0 0 0 {L*MM} 0 0" size="0.003"/>
      <site name="OA_end" pos="{L*MM} 0 0" size="0.0035"/>
    </body>
  </worldbody>
</mujoco>"""
m = mujoco.MjModel.from_xml_string(xml_one)
d = mujoco.MjData(m)
d.qpos[0] = np.radians(30)              # rotate the link 30 deg
mujoco.mj_forward(m, d)                 # compute positions from qpos
print("single link end at 30 deg (mm):", d.site("OA_end").xpos / MM)



def link(name, length_mm, rgba, children=""):
    l = length_mm * MM
    return f"""
      <joint name="{name}" type="hinge" axis="0 0 1"/>
      <geom type="capsule" fromto="0 0 0 {l} 0 0" size="0.003" rgba="{rgba}"/>
      <site name="{name}_end" pos="{l} 0 0" size="0.0035" rgba="1 1 1 1"/>
      {children}"""


blue, orange, green = "0.12 0.47 0.71 1", "1 0.5 0.05 1", "0.17 0.63 0.17 1"

tree = f"""
    <body name="OA" pos="0 0 0">{link("OA", L, blue, f'''
      <body name="DA" pos="{L*MM} 0 0">{link("DA", s, green)}</body>''')}
    </body>

    <body name="OC" pos="0 0 0">{link("OC", L, blue, f'''
      <body name="CD" pos="{L*MM} 0 0">{link("CD", s, green)}</body>''')}
    </body>

    <body name="EB" pos="{r*MM} 0 0">{link("EB", r, orange, f'''
      <body name="AB" pos="{r*MM} 0 0">{link("AB", s, green)}</body>
      <body name="BC" pos="{r*MM} 0 0">{link("BC", s, green)}</body>''')}
    </body>"""

m = mujoco.MjModel.from_xml_string(f"<mujoco><worldbody>{tree}</worldbody></mujoco>")
print("bodies:", m.nbody, " joints:", m.njnt)



equality = """
  <equality>
    <connect name="cut_A" site1="AB_end" site2="OA_end"/>
    <connect name="cut_C" site1="BC_end" site2="OC_end"/>
    <connect name="cut_D" site1="DA_end" site2="CD_end"/>
  </equality>"""

CUTS = [("AB_end", "OA_end"), ("BC_end", "OC_end"), ("DA_end", "CD_end")]


def gaps(d):
    return [np.linalg.norm(d.site(a).xpos - d.site(b).xpos) / MM for a, b in CUTS]


m = mujoco.MjModel.from_xml_string(f"<mujoco><worldbody>{tree}</worldbody>{equality}</mujoco>")
d = mujoco.MjData(m)
mujoco.mj_forward(m, d)
print("equality constraints:", m.neq)
print("gaps at qpos = 0 (mm):", np.round(gaps(d), 2))



x_line = (L**2 - s**2) / (2 * r) # predicted straight line, for the video overlay


def build_xml():
    return f"""
<mujoco model="peaucellier_lipkin">
  <option gravity="0 0 0" timestep="0.0002" integrator="implicitfast">
    <flag contact="disable"/>
  </option>
  <default>
    <joint damping="0.002" armature="1e-6"/>
    <geom density="1000" contype="0" conaffinity="0"/>
    <equality solref="0.001 1" solimp="0.99 0.999 0.0001"/>
  </default>
  <visual><global offwidth="960" offheight="720"/></visual>
  <worldbody>
    <light pos="0.08 0 1" dir="0 0 -1"/>
    <camera name="top" pos="0.075 0 0.45" xyaxes="1 0 0 0 1 0" fovy="45"/>
    <geom type="plane" size="0.3 0.3 0.01" pos="0.075 0 -0.01" rgba="0.95 0.95 0.95 1"/>
    <geom type="box" size="0.0006 0.15 0.0005" pos="{x_line*MM} 0 -0.004" rgba="0.85 0.1 0.1 1"/>
    {tree}
  </worldbody>
  {equality}
  <actuator>
    <position name="crank" joint="EB" kp="20" kv="0.02"/>
  </actuator>
</mujoco>"""


m = mujoco.MjModel.from_xml_string(build_xml())
d = mujoco.MjData(m)
print("actuators:", m.nu, " timestep:", m.opt.timestep)



def qadr(name):
    return m.joint(name).qposadr[0]


def tracer(d):
    return 0.5 * (d.site("DA_end").xpos + d.site("CD_end").xpos)[:2] / MM


# These only choose the branch: A above OB, C below, rhombus open.
guess = dict(OA=20, OC=-20, EB=0, AB=70, BC=-70, DA=-90, CD=90)
for name, deg in guess.items():
    d.qpos[qadr(name)] = np.radians(deg)
mujoco.mj_forward(m, d)
print("gaps before (mm):", np.round(gaps(d), 2))

d.ctrl[0] = 0.0 # hold the crank at theta = 0
for _ in range(int(1.0 / m.opt.timestep)): # simulate 1 s
    mujoco.mj_step(m, d)
print("gaps after  (mm):", ["%.1e" % g for g in gaps(d)])
print("D at theta=0 (mm):", tracer(d))



dt = m.opt.timestep
frames, path = [], []
renderer = mujoco.Renderer(m, 720, 960) # for the video


def step(n):
    for _ in range(n):
        mujoco.mj_step(m, d)
        if int(round(d.time / dt)) % 167 == 0: # about 30 fps
            renderer.update_scene(d, camera="top")
            frames.append(renderer.render())
            path.append([d.time, *tracer(d), *gaps(d)])


def go_to(target_deg, speed_dps=40, hold_s=0.4):
    start, target = d.ctrl[0], np.radians(target_deg)
    n = max(1, int(abs(target - start) / np.radians(speed_dps) / dt))
    for k in range(1, n + 1):
        d.ctrl[0] = start + (target - start) * k / n  # slow ramp
        step(1)
    step(int(hold_s / dt)) # settle


TEST = [-80, -60, -40, 0, 40, 60, 80] # test angles (deg)
go_to(TEST[0], hold_s=0.6)
log = []
for t in TEST:
    go_to(t)
    th = np.degrees(d.qpos[qadr("EB")])  # measured angle
    Dx, Dy = tracer(d)
    log.append([th, Dx, Dy, *gaps(d), np.abs(d.qvel).max()])
    print(f"theta {th:8.3f}  D = ({Dx:9.4f}, {Dy:9.4f})  max gap {max(gaps(d)):.1e}  "
          f"|qvel| {np.abs(d.qvel).max():.0e}")
log = np.array(log)



path = np.array(path)
moving = path[:, 0] > 1.0
print("max gap while MOVING (mm):", path[moving, 3:6].max())
print("max gap when SETTLED (mm):", log[:, 3:6].max())
np.savetxt("d8_mujoco_log.csv", log, delimiter=",",
           header="theta_deg,Dx_mm,Dy_mm,gapA,gapC,gapD,max_qvel")



print(f"{len(frames)} frames = {len(frames)/30:.1f} s")
imageio.mimsave("d8_video.mp4", frames, fps=30, macro_block_size=8)
print("saved d8_video.mp4")
