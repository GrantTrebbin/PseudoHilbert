# PseudoHilbert

This is a python module to create Pseudo Hilbert curves to cover arbitrary rectanglular regions. It is an implementation of the paper

"A Pseudo-Hilbert Scan for Arbitrarily-Sized Arrays" by Jian ZHANG, Sei-ichiro KAMATO, Yoshifumi UESHIGE

Improvements from later work by one of the Authors are also included.

The PseudoHilbert object will have two attributes for indexing.
coordinate_to_index Allows the lookup of the index when supplied with x and y coordinates
index_to_coordinate Allows the lookup of the coordinate when supllied with the index

It's important to note that the curve will start in a corner but won't necessarily end in one.
It will however end near a corner.  Every cell will be covered though.

It works and is reasonalbly fast, but I can already think of better ways to implement it.
