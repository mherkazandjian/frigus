Solve your own system
=====================

This page is intended to guide a user to define their own system of levels
in the framework of Frigus, such that the following would be possible:

.. code-block:: python

     from frigus.readers import DataLoader
     species = DataLoader().load('my_custom_species')
     cooling_rate = cooling_rate_at_steady_state(species, T_kin, T_rad, nc)
