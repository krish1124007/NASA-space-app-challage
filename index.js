// Animation loop
function animate() {
    requestAnimationFrame(animate);

    // Only animate if data is loaded
    if (!isDataLoaded || !asteroidData) {
        renderer.render(scene, camera);
        return;
    }

    if (!isPaused) {
        // Update Earth position using the same orbital calculation as the orbit line
        const earthM = (time / 365.25) * 360; // Mean anomaly in degrees
        const earthPos = orbitalToCartesian(
            earthOrbitRadius,  // semi-major axis
            earthEccentricity, // eccentricity
            0,                 // inclination (Earth's reference plane)
            0,                 // RAAN
            0,                 // argument of periapsis
            earthM             // mean anomaly
        );
        
        // Move the container (orbital position)
        earthContainer.position.x = earthPos.x;
        earthContainer.position.y = earthPos.y;
        earthContainer.position.z = earthPos.z;
        
        // Rotate Earth on its tilted axis (daily rotation)
        earth.rotation.y += 0.01;

        // Update asteroid
        const period = Math.pow(asteroidData.a, 1.5) * 365.25;
        const M = asteroidData.M0 + (time / period) * 360;
        const scaledA = asteroidData.a * AU_TO_UNITS;
        const astPos = orbitalToCartesian(scaledA, asteroidData.e, asteroidData.i,
            asteroidData.omega, asteroidData.w, M);
        
        asteroid.position.set(astPos.x, astPos.y, astPos.z);
        asteroid.rotation.y += 0.02;

        // Check for collision (use earthContainer position)
        const distToEarth = Math.sqrt(
            Math.pow(astPos.x - earthContainer.position.x, 2) +
            Math.pow(astPos.y - earthContainer.position.y, 2) +
            Math.pow(astPos.z - earthContainer.position.z, 2)
        );
        
        const distToEarthKm = (distToEarth / AU_TO_UNITS) * AU_TO_KM;
        
        if (distToEarthKm <= EARTH_RADIUS_KM && collisionDetected) {
            isPaused = true;
            document.getElementById('collision-impact').style.display = 'block';
        }

        // Update HUD
        const velocity = calculateOrbitalVelocity(astPos.r / AU_TO_UNITS, asteroidData.a);
        
        document.getElementById('hud-time').textContent = `Day ${Math.floor(time)}`;
        document.getElementById('hud-anomaly').textContent = (astPos.nu * 180 / Math.PI).toFixed(2) + '°';
        document.getElementById('hud-speed').textContent = velocity.toFixed(2) + ' km/s';
        document.getElementById('hud-distance').textContent = ((distToEarth / AU_TO_UNITS) * AU_TO_KM / 1e6).toFixed(2) + ' million km';
        document.getElementById('hud-radius').textContent = (astPos.r / AU_TO_UNITS).toFixed(3) + ' AU';

        sun.rotation.y += 0.005;
        time += timeSpeed;
    }

    renderer.render(scene, camera);
}