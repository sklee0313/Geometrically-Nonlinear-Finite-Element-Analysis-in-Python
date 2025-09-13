# Geometrically Nonlinear Finite Element Analysis in Python

Available examples are:  

**1. Slender Beam Bending Example**  
   A two-dimensional slender beam is subject to a point force that creates bending-type large displacements.  
   The 8-node quadrilateral element is used to effectively capture the bending, and the Total Lagrangian formulation is employed with plane stress conditions and the Saint Venant–Kirchhoff model.  
   The problem is solved in [here](https://doi.org/10.1016/0045-7949(77)90027-X) and [here](https://doi.org/10.1108/EC-10-2019-0478).

**2. Slender Beam Bending UL Example**  
   The same problem as in the Slender Beam Bending Example is solved but using the Updated Lagrangian formulation. Small strains are assumed, so no transformation is applied to the material model.

**3. Circular Arc Beam Example**  
   A three-dimensional arc beam is subject to a tip force that creates bending and torsion.  
   The domain is discreized using the 20-node brick element, and the Total Lagrangian formulation is employed with the Saint Venant–Kirchhoff model.
   The problem was introduced by [this](https://onlinelibrary.wiley.com/doi/10.1002/nme.1620140703).
