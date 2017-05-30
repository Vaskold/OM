using System;

namespace GoingDownForReal
{
    internal class Program
    {
        public struct Vector
        {
            public double X1;
            public double X2;

            public static Vector operator -(Vector x, Vector y)
            {
                Vector z;
                z.X1 = x.X1 - y.X1;
                z.X2 = x.X2 - y.X2;
                return z;
            }

            public static Vector operator +(Vector x, Vector y)
            {
                Vector z;
                z.X1 = x.X1 + y.X1;
                z.X2 = x.X2 + y.X2;
                return z;
            }

            public static Vector operator *(double a, Vector y)
            {
                Vector z;
                z.X1 = a * y.X1;
                z.X2 = a * y.X2;
                return z;
            }
        }

        private static Vector Mul(double[][] a, Vector b)
        {
            Vector r;
            r.X1 = a[0][0] * b.X1 + a[0][1] * b.X2;
            r.X2 = a[1][0] * b.X1 + a[1][1] * b.X2;
            return r;
        }

        private static void Fill(double[][] a, double a1, double a2, double a3, double a4)
        {
            a[0][0] = a1;
            a[0][1] = a2;
            a[1][0] = a3;
            a[1][1] = a4;
        }

        private static double F(Vector x) // your f(x)
        {
            return
                x.X1 * x.X1 +
                7 * x.X2 * x.X2 +
                x.X1 * x.X2 +
                x.X1;
        }

        private static double Ff(Vector xk, Vector grad, double x) // f(a)
        {
            return
                Math.Pow(xk.X1 - x * grad.X1, 2) +
                7 * Math.Pow(xk.X2 - x * grad.X2, 2) +
                (xk.X1 - x * grad.X1) * (xk.X2 - x * grad.X2) +
                (xk.X1 - x * grad.X1);

            /*
       return 
            6 * (-1) * grad.x1 * (xk.x1 - x * grad.x1) + 
            2 * (-1) * grad.x2 * (xk.x2 - x * grad.x2) - 
            ((-1) * grad.x1 * (xk.x2 - x * grad.x2) + (-1) * grad.x2 * (xk.x1 - x * grad.x1)) - 
            4 * (-1) * grad.x1;
            */
        }

        private static Vector Grad(Vector xk) // return grad(f(xk))
        {
            Vector grd;
            grd.X1 = 2 * xk.X1 + xk.X2 + 1;
            grd.X2 = 14 * xk.X2 + xk.X1;
            return grd;
        }

        public static double Gs(Vector xk, Vector grad)
        {
            double a = 0, b = 1;
            const double del = 0.0001;
            const double eps = 0.002;

            while (Math.Abs(b - a) > eps)
            {
                var x = (a + b) / 2 - del;
                var y = (a + b) / 2 + del;

                if (Ff(xk, grad, x) < Ff(xk, grad, y))
                {
                    b = y;
                }
                else if (Ff(xk, grad, x) > Ff(xk, grad, y))
                {
                    a = x;
                }
                else
                {
                    a = x;
                    b = y;
                }
            }
            return (a + b) / 2;
        } // return a

        private static double Norma(Vector x)
        {
            return Math.Sqrt(Math.Pow(x.X1, 2) + Math.Pow(x.X2, 2));
        }

        private static double Sqr(Vector x)
        {
            return Math.Pow(x.X1, 2) + Math.Pow(x.X2, 2);
        }

        private static void Main()
        {
            var a = new double[2][];
            for (var j = 0; j < 2; j++)
                a[j] = new double[2];

            const double e1 = 0.01;
            const double e2 = 0.01;

            int i = 1, cond = 0;
            //bool cond = true;

            Vector x0;
            x0.X1 = -1.1;
            x0.X2 = 1.1;

            var f0 = F(x0);
            var f1 = f0;
            var p = -1 * Grad(x0);
            var gradSqr = Sqr(p);

            Fill(a, 0.519, -0.037, -0.037, 0.074);

            while (Norma(Grad(x0)) > e1) // && cond < 2)
            {
                Console.WriteLine("\nStep {0}\n", i);
                ++i;

                Console.WriteLine("x{4} = ({0}; {1}), grad = ({2} {3})",
                    x0.X1, x0.X2, Grad(x0).X1, Grad(x0).X2, i - 2);

                Console.WriteLine("|Grad(f(x{0}))| = {1}", i - 2, Norma(Grad(x0)));

                Console.WriteLine("f(x{0}) = {1}", i - 2, F(x0));

                //Console.WriteLine("a{1} = {0}", Gs(x0, Grad(x0)), i - 2); show alpha

                //var xk = x0 - Gs(x0, Grad(x0)) * Grad(x0);

                var xk = x0 + Gs(x0, Grad(x0)) * p;
                var newGrad = -1 * Grad(xk);
                var newGradScr = Sqr(newGrad);
                var b = newGradScr / gradSqr;
                p = newGrad + b * p;

                //var xk = x0 - Mul(a, Grad(x0));

                Console.WriteLine("f(x{0}) = {1}", i - 1, F(xk));

                Console.WriteLine("x{2} = ({0}; {1})", xk.X1, xk.X2, i - 1);

                if (F(xk) - F(x0) < e2 && Norma(xk - x0) < e2)
                {
                    ++cond;
                    //Console.WriteLine("{0} {1}", F(xk) - F(x0), Norma(xk - x0));
                    Console.WriteLine("{0} check is completed", cond);
                    if (Norma(Grad(xk)) < e1)
                        Console.WriteLine("|Grad(f(x{0}))| = {1}, complete", i - 1, Norma(Grad(xk)));
                }

                x0 = xk;
                f0 = f1;
                f1 = F(xk);
                gradSqr = newGradScr;
            }

            Console.ReadKey();
        }
    }
}