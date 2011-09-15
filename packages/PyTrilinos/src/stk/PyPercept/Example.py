import sys
sys.path.append("/scratch/srkenno/Trilinos-BUILDS/build11-090711/packages/PyTrilinos/src/stk/PyPercept")

from mpi4py import MPI
from PerceptMesh import *
import unittest
from numpy import *
import subprocess

class UseCases(unittest.TestCase):

    def test_use_case_2(self):

        pMesh = PerceptMesh()                          # create an empty PerceptMesh
        pMesh.open("exodus_files/tet-mesh.e")          # open the mesh, but don't commit its meta data

        uniform_refiner = Refiner(pMesh, TET4_TET4_8)  # define a Refiner on the mesh
        pMesh.commit()                                 # commit the mesh

        i = 0
        while i < 3:
            uniform_refiner.doBreak()                  # refine the mesh 3 times
            i = i + 1

        pMesh.saveAs("tet-mesh-refined-3-times.e")     # save in exodus format

    def test_fieldFunction_demo_2(self):   
        eMesh = PerceptMesh()
        eMesh.newMesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"))      # use a fixture to generate a 3x3x3 hex mesh

        vectorDimension = 0
        # add a field
        eMesh.addField("coords_mag_field", FEMMetaData.NODE_RANK, vectorDimension)  
        eMesh.commit()

        f_coords = eMesh.getField("coordinates")                # get pre-existing field
        coords_mag_field = eMesh.getField("coords_mag_field")   # get the field we just created

        ff_coords = FieldFunction("ff_coords", f_coords, eMesh, 3, 3)  # define a field function
        evalVec3Print(0.1,0.1,0.1,0.0,ff_coords)  # evaluate and print the field function a point {0.1, 0.1, 0.1} time=0.0

        coords_mag_sf = StringFunction("sqrt(x*x + y*y + z*z)" , "coords_mag_sf", 3, 1)  # define coordinate magnitude function
        x = 0.123
        y = 0.234
        z = 0.345
        vv = sqrt(x*x + y*y + z*z)
        v1 = evalFunc(x,y,z,0,coords_mag_sf)
        print "vv = ", vv, "== v1 = ", v1
        self.assertEqual(vv, v1)               # ensure correctness of string function

        # define a field function
        coords_mag_field_function = FieldFunction("coords_mag_field_function", coords_mag_field, eMesh, 3, 1)

        # interpolate the function onto the mesh
        coords_mag_field_function.interpolateFrom(coords_mag_sf)

        eMesh.saveAs("./cubehex8_withCoordMag_out.e")

        # demonstrate how to usa an alias
        ff_coords.addAlias("mc")

        sfcm = StringFunction("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", "sfcm", 3, 1)

    def test_wedge6_wedge18_enrich(self):
        pm = MPI.COMM_WORLD
        p_size = parallel_machine_size(pm)
        if p_size == 1:
            wedgeFixture = WedgeFixture()
            bulk = wedgeFixture.createMesh(MPI.COMM_WORLD, 4,3,2,0,1,0,1,0,1,"")   # create stk::mesh::BulkData from wedge fixture
            eMesh = PerceptMesh(wedgeFixture.getMetaData(), bulk, False)           # adopt bulk data
            scalarDimension = 0
            proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension)
            breaker = Refiner(eMesh, WEDGE6_WEDGE15_1, proc_rank_field)
            eMesh.commit()
            wedgeFixture.createBulkAfterMetaCommit(MPI.COMM_WORLD)         # generate the mesh
            breaker.doBreak()                                              # refine
            eMesh.saveAs("./wedge6-15.e")                                  # save


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(UseCases)
    unittest.TextTestRunner(verbosity=2).run(suite)
