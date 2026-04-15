"""
Struct.py

Written by: Steven Chen

Converts from dictionary to object.

CHANGELOG:

    01/21/2014 - Converts a dictionary to a object-like format (struct)
    
"""

class Struct:
    
    def __init__(self, **entries):
        for key, value in entries.items():
            try:
                if isinstance(value, dict):    
                    value2 = Struct(**value)
                else:
                    value2 = value
            except:
                value2 = value
            self.__dict__[key] = value2