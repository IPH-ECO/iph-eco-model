# iph-eco-model
IPH-ECO source code
IPH-ECO is a three-dimensional complex dynamic model of aquatic ecosystems such as lakes,
reservoirs and estuaries. The model consists of several partial differential equations being
divided in three modules: (a) a detailed hydrodynamic module, describing quantitative flows
and water level; (b) a nutrient module, which deals with nutrient transport mechanisms;
and (c) a biological module, which describes whole aquatic food-web interactions.
The differential equations are solved numerically by applying an efficient semi-implicit finite differences
method in a three-dimensional regular grid.
For an extended model like IPH-ECO it is important think about how to control the complexity with the numerous parameters.
An effective way to control complexity is to implement flexibility in the model design in such way that it is possible to
switch between different modes of complexity. Therefore IPH-ECO includes a graphical user-friendly interface (GUI)
for MS Windows environment with a flexible design to vary the complexity of the model. The GUI was developed in Visual Basic
which is an event-driven programming language with an integrated development environment (IDE) for its Component Object Models (COM),
making it easier the building of interfaces for applications in any OS. As the simulation model was implemented in Visual FORTRAN,
which is connected with GUI by a FORTRAN Dynamic Link Library (DLL).