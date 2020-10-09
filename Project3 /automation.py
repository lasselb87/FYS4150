import os

# Compile programs and run tests
os.system("make all")

# Earth-Sun system
os.system("./twoBody.x 5 10000 10")
os.system("./benchmark.x")

# Escape velocity of Earth
os.system("./escapeVelocity.x")

# Earth-Jupiter-Sun system
os.system("./threeBody.x")

# Orbits of all planets in the Solar System
os.system("./allPlanets.x")

# Perihelion precession of Mercury
os.system("./relativistic.x")

# Clean
os.system("make clean")
