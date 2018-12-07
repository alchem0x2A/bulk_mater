import numpy
import os
import matplotlib.pyplot as plt
import unittest

class Test(unittest.TestCase):
    def test_eps(self):
        base_dir = os.path.dirname(__file__)
        data = numpy.load(os.path.join(base_dir,
                                       "../../tmp/Si-eps-test/eps_rpa.npz"))
        freq = data["frequencies"]
        epsx = data["eps_x"]
        epsy = data["eps_y"]
        epsz = data["eps_z"]
        self.assertLess(numpy.abs(numpy.mean(epsx - epsy)), 1e-4)
        self.assertLess(numpy.abs(numpy.mean(epsx - epsz)), 1e-4)
        plt.plot(freq, epsx.real, label="epsilon_1")
        plt.plot(freq, epsx.imag, label="epsilon_2")
        plt.xlabel("E (eV)")
        plt.ylabel("Epsilon")
        # plt.xlim(0, 10)
        plt.legend()
        plt.show()
        # plt.savefig("dielectric-test.png")
        return

if __name__ == "__main__":
    unittest.main()
