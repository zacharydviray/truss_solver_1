Terminology:

NoN: Number of Nodes
NPE: Nodes per Element
NoE: Number of Elements
PD: Problem Dimension
EL: Element List
ENL: Extended Node List (Size of NoN x 6*PD)

NL (NoN x PD) -> These refer to the coordinates of the nodes themselves

EL (Number of Elements x Nodes per Element) -> These refer to what nodes each element will consist of

*This assumes that each element will consist of two nodes, but that's chill. 
DorN (D/N): Matrix that represents assigned boundary conditions (-1 if fixed, 1 if free)

ENL will consist of nl, DorN, Local Degrees of Freedom, Global Degrees of Freedom, Displacement, and Force

nl is a mini node list that represent the node number of the selected element
