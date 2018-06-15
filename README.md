# Graphics Engine
## by Ricky Chen - Period 10

### Implemented Commands
* LIGHT -- modified from specification to support multiple lights (a default light is supplied if none is defined).
* * light <name> <location x, y, z> <color r, g, b>
* AMBIENT
* * ambient <color r, g, b>
* SHADING -- supports flat, gouraud, phong
* * shading <flat, gouraud, phong>

## Premade Animation
Multiple lights and phong shading can be seen in action by doing

$ make

$ ./mdl test_ball.mdl

$ display ball.gif
