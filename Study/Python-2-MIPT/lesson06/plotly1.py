import plotly.graph_objs as go
import numpy as np

def cubes(size, pos_x, pos_y, pos_z, color):
    # create points
    x, y, z = np.meshgrid(
        np.linspace(pos_x-size/2, pos_x+size/2, 2), 
        np.linspace(pos_y-size/2, pos_y+size/2, 2), 
        np.linspace(pos_z-size/2, pos_z+size/2, 2),
    )
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()
    
    return go.Mesh3d(x=x, y=y, z=z, alphahull=1, flatshading=True, color=color, lighting={'diffuse': 0.1, 'specular': 2.0, 'roughness': 0.5})

fig = go.Figure()

bigsize = 20
size = 5

# add outer cube
fig.add_trace(cubes(bigsize, 0,0,0, 'rgba(255,100,0,0.1)'))

# add inner center cube
fig.add_trace(cubes(size, 0,0,0, 'rgba(100,0,100,0.1)'))

# add inner cubes
fig.add_trace(cubes(size, size,0,0, 'rgba(100,0,100,0.1)'))
fig.show()
