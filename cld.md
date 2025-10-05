1️⃣ Asteroid Orbit aur Trajectory

Har asteroid ka orbit hota hai jo sun ke around hota hai. Orbit ke parameters hote hain:

Semi-major axis (a) – orbit ka size

Eccentricity (e) – orbit kitni elliptical hai

Inclination (i) – orbit ka tilt Earth ke plane ke respect

Longitude of ascending node (Ω) – orbit ka orientation

Argument of periapsis (ω) – closest point Sun se

Mean anomaly (M) – asteroid ki current position orbit me

Ye data NASA APIs (JPL Small-Body Database) se milta hai.

2️⃣ Earth ke position ke sath compare karna

Earth ka bhi orbit hota hai sun ke around. Ab hume asteroid ki trajectory aur Earth ke orbit ko time ke function ke roop me plot karna hota hai:

Asteroid ke orbit ko ke 3D Cartesian coordinates (x, y, z) me convert karte hain using orbital elements.

Earth ke orbit ke liye bhi same karte hain.

Formulas (approximate):
For asteroid position at time t:

𝑥
=
𝑟
(
cos
⁡
(
Ω
)
cos
⁡
(
𝜃
+
𝜔
)
−
sin
⁡
(
Ω
)
sin
⁡
(
𝜃
+
𝜔
)
cos
⁡
(
𝑖
)
)
x=r(cos(Ω)cos(θ+ω)−sin(Ω)sin(θ+ω)cos(i))
𝑦
=
𝑟
(
sin
⁡
(
Ω
)
cos
⁡
(
𝜃
+
𝜔
)
+
cos
⁡
(
Ω
)
sin
⁡
(
𝜃
+
𝜔
)
cos
⁡
(
𝑖
)
)
y=r(sin(Ω)cos(θ+ω)+cos(Ω)sin(θ+ω)cos(i))
𝑧
=
𝑟
(
sin
⁡
(
𝜃
+
𝜔
)
sin
⁡
(
𝑖
)
)
z=r(sin(θ+ω)sin(i))

𝑟
r = distance from Sun at that point of orbit

𝜃
θ = true anomaly at time t

3️⃣ Miss distance calculate karna

Miss distance = Earth aur asteroid ke closest approach ka distance.

Calculate asteroid ka x, y, z at time t

Calculate Earth ka x, y, z at same time t

Distance formula use karo:

𝑑
=
(
𝑥
𝑎
−
𝑥
𝑒
)
2
+
(
𝑦
𝑎
−
𝑦
𝑒
)
2
+
(
𝑧
𝑎
−
𝑧
𝑒
)
2
d=
(x
a
	​

−x
e
	​

)
2
+(y
a
	​

−y
e
	​

)
2
+(z
a
	​

−z
e
	​

)
2
	​


Agar 
𝑑
≤
𝑅
𝐸
𝑎
𝑟
𝑡
ℎ
+
𝑠
𝑎
𝑓
𝑒
𝑡
𝑦
_
𝑚
𝑎
𝑟
𝑔
𝑖
𝑛
d≤R
Earth
	​

+safety_margin, matlab asteroid Earth se takrayega

Agar 
𝑑
>
𝑅
𝐸
𝑎
𝑟
𝑡
ℎ
d>R
Earth
	​

, matlab Earth ke paas se hoga par takrayega nahi

NASA aur ESA scientists ye calculate karte hain orbit propagation ke algorithms se, jaise Runge-Kutta numerical integration, because orbit ke forces complex hote hain (Sun, planets, radiation, etc).

4️⃣ Hazard assessment

Agar asteroid >140 m aur miss distance <0.05 AU (~7.5 million km), toh ye Potentially Hazardous Asteroid (PHA) hota hai.

🔹 Example

Suppose asteroid ka semi-major axis 1.2 AU, eccentricity 0.1

Inclination 5°

Hum calculate karte hain asteroid aur Earth ka distance 2026-01-01 ko:

Earth x,y,z = (0.9, 0.2, 0) AU

Asteroid x,y,z = (0.95, 0.22, 0.01) AU

Distance 
𝑑
=
𝑠
𝑞
𝑟
𝑡
(
(
0.95
−
0.9
)
2
+
(
0.22
−
0.2
)
2
+
(
0.01
−
0
)
2
)
≈
0.057
𝐴
𝑈
d=sqrt((0.95−0.9)
2
+(0.22−0.2)
2
+(0.01−0)
2
)≈0.057AU → safe, no collision

✅ Summary (Simple English)

Get asteroid orbital elements (a, e, i, Ω, ω, M)

Convert orbit into 3D Cartesian coordinates over time

Calculate Earth's position at the same time

Compute distance between asteroid and Earth

If distance < Earth radius + safety margin → collision, else safe