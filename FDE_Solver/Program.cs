using System;
using System.Threading;
using System.Threading.Tasks;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace FDE_Solver
{
	class MainClass
	{
		public delegate double ForcingFunction(double x, double y);


		public static void Main (string[] args)
		{
			Console.WriteLine ("Hello World!");
			Console.WriteLine (SpecialFunctions.Factorial(5));
			double[] y = Compute (0.5, new double[] { 0 }, Math.PI, 1000, 2, new ForcingFunction (ff));
			StreamWriter sw = new StreamWriter ("data.csv");
			foreach (double v in y) {
				sw.WriteLine (v.ToString ());
			}
			sw.Close ();

		}
		public static double ff(double x, double y)
		{
			return Math.Tan (x);
		}


		public static double a(int mu, double alpha)
		{
			return (Math.Pow(mu + 2, alpha + 1) - 2 * Math.Pow(mu + 1, alpha + 1) + Math.Pow(mu, alpha + 1))/SpecialFunctions.Gamma(alpha + 2);
		}

		public static double c(int mu,  double alpha)
		{
			return (Math.Pow(mu, alpha+1) - (mu - alpha)*Math.Pow(mu + 1, alpha)) / SpecialFunctions.Gamma(alpha + 2);
		}

		public static double b(int mu, double alpha)
		{
			return (Math.Pow (mu + 1, alpha) - Math.Pow (mu, alpha)) / SpecialFunctions.Gamma (alpha + 1);
		}

		public static double I_1(int j, double alpha, double[] y_0_diffs, double[] x)
		{
			double value = 0;
			for (int k = 0; k <= Math.Ceiling (alpha) - 1; k++) {
				value += Math.Pow (x [j + 1], k) / SpecialFunctions.Factorial (k) * y_0_diffs [k];
			}
			return value;
		}

		public static double H_p(int j, int ell, int p, double[] x, double[] y, double alpha, ForcingFunction f)
		{
			double value = 0;
			for (int k = 0; k <= (ell - 1) * p; k++) {
				value += b (j - k, alpha) * f (x [k], y [k]);
			}
			return value;
		}

		public static double L_p(int j, int ell, int p, double[] x, double[] y, double alpha, ForcingFunction f) 
		{
			double value = 0;
			for (int k = (ell - 1) * p + 1; k <= j; k++) {
				value += b (j - k, alpha) * f (x [k], y [k]);
			}
			return value;
		}

		public static double H(int j, int ell, int p, double[] x, double[] y, double alpha, ForcingFunction f)
		{
			double value = 0;
			value += c (j, alpha) + f (x [0], y [0]);
			for (int k = 1; k <= (ell - 1) * p; k++) {
				value += a (j - k, alpha) * f (x [k], y [k]);
			}
			return value;
		}

		public static double L(int j, int ell, int p, double[] x, double[] y, double alpha, ForcingFunction f, double y_p_1)
		{
			double value = 0;
			for (int k = (ell - 1)*p + 1; k <= j; k++) {
				value += a(j-k, alpha) * f(x[k], y[k]);
			}
			value += f(x[j+1], y_p_1) / SpecialFunctions.Gamma(alpha + 2);
			return value;
		}

		public static double[] Compute(double alpha, double[] y_0_diffs, double T, int N, int p, ForcingFunction f)
		{
			double[] x = new double[N];
			double[] y = new double[N];
			double[] y_p = new double[N];
			y [0] = y_0_diffs [0];
			double h = T / N;
			for (int i = 0; i < N; i++)
			{
				x [i] = h * i;
			}
			for (int ell = 1; ell <= Math.Ceiling ((double)N / (double)p); ell++) {
				Task<double> taskSum_p = null;
				Task<double> taskSum = null;
				for (int i = 0; i < p && ((ell - 1) * p) + i < N - 1; i++) {
					int j = ((ell - 1) * p) + i;
					Task<double> taskI = Task.Factory.StartNew (() => I_1 (j, alpha, y_0_diffs, x));

					Task<double> taskH_p = Task.Factory.StartNew (() => H_p (j, ell, p, x, y, alpha, f));
					Task<double> taskH = Task.Factory.StartNew (() => H (j, ell, p, x, y, alpha, f));
					Task<double> taskL_p = null;
					if (taskSum != null) {
						taskL_p = taskSum.ContinueWith ((t) => L_p (j, ell, p, x, y, alpha, f));
					} else {
						taskL_p = Task.Factory.StartNew (() => L_p (j, ell, p, x, y, alpha, f));
					}
					taskSum_p = Task.Factory.ContinueWhenAll(new [] { taskL_p, taskH_p, taskI }, (ts) => y_p[j + 1] = taskI.Result + Math.Pow(h, alpha) * ( taskH_p.Result + taskL_p.Result ) );
					Task<double> taskL = taskSum_p.ContinueWith ((t) => L (j, ell, p, x, y, alpha, f, y_p [j + 1]));
					taskSum = Task.Factory.ContinueWhenAll(new [] { taskH, taskL, taskI }, (ts) => y[j+1] = taskI.Result + Math.Pow(h, alpha) * (taskH.Result + taskL.Result ));
				}
				Console.WriteLine (ell.ToString () + " / " + Math.Ceiling ((double)N / (double)p).ToString ());
				taskSum.Wait ();
			}
			return y;
		}



	}
}
