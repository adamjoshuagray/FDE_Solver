using System;

namespace FDE_Solver
{
	public class SpecialFunctions
	{
		public static double Gamma(double z)
		{
			double g = 7;
			double[] p =  new double[] { 0.99999999999980993, 676.5203681218851, -1259.1392167224028,
										 771.32342877765313, -176.61502916214059, 12.507343278686905,
				-0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };
			if (z < 0.5) {
				return Math.PI / (Math.Sin (Math.PI * z) * Gamma (1 - z));
			} else {
				z -= 1;
				double x = p [0];
				for (int i = 1; i < g + 2; i++)
				{
					x += p [i] / (z + i);
				}
				double t = z + g + 0.5;
				return Math.Sqrt (2 * Math.PI) * Math.Pow (t, z + 0.5) * Math.Exp (-t) * x;
			}
		}
		public static int Factorial(int k)
		{
			int value = 1;
			for (int i = 1; i <= k; i++) {
				value *= i;
			}
			return value;
		}
	}
}

