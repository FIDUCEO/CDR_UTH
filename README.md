# CDR_UTH

Development code for UTH CDR

    Copyright (C) 2019 Theresa Lang University of Hamburg
    This code was developed for the EC project “Fidelity and Uncertainty in
    Climate Data Records from Earth Observations (FIDUCEO)”.
    Grant Agreement: 638822
    
    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your option)
    any later version.
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.
    A copy of the GNU General Public License should have been supplied along
    with this program; if not, see http://www.gnu.org/licenses/


The folder processing contains the Python code to process the 1c FIDUCEO Microwave FCDR to the l3 UTH CDR.

Note that you will need additional code to execute the processing: The CDR-generator code makes use of the open source code "typhon" (available through http://www.radiativetransfer.org/tools/) and the FIDUCEO NetCDF writer (available through https://github.com/FIDUCEO/FCDRTools/tree/master/fiduceo/cdr)

An overview of the code and its usage is provided in overview_CDR_processing.txt in processing.

A description of the CDR processing chain can be found in the UTH CDR Product User Guide V1.0.
