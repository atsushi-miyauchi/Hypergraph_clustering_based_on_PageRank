using System;

using MathNet.Numerics.LinearAlgebra;
using System.Collections.Generic;
using System.Linq;

namespace SubmodularHeatEquation
{
    class MainClass
    {
        public static double Sqr(double v) { return v * v; }
        const double dt = 0.01;


        public static Vector<double> Integrate(Func<double, Vector<double>> f, double T)
        {
            var dummy = f(0.0);
            var res = CreateVector.Dense<double>(dummy.Count);

            var t = 0.0;
            while (t < T)
            {
                double nt = Math.Min(T - t, dt);
                var v = f(t);
                res += v * nt;
                t += nt;
            }
            return res;
        }


        // for Proposed_local
        public static void Proposed_local(string fn, int v_init)
        {
            var H = Hypergraph.Open(fn);

            var time = new System.Diagnostics.Stopwatch();
            time.Start();

            int n = H.n;
            int m = H.m;

            const double eps = 0.9;

            const double dt = 1.0;
            const double T = 30.0;

            var A_cand = new List<double>();
            for (int i = 0; i <= Math.Log(n * m) / Math.Log(1 + eps); i++)
            {
                A_cand.Add(Math.Pow(1 + eps, i) / (n * m));
            }

            var edge_size = new Dictionary<int, int>();
            for (int eid = 0; eid < H.m; eid++)
            {
                edge_size.Add(eid, H.ID_rev[eid].Count());
            }

            double min_conductance = double.MaxValue;

            foreach (double alpha in A_cand)
            {
                var vec = CreateVector.Dense<double>(n);

                vec[v_init] = 1.0;
                
                vec = Hypergraph.Simulate(H, vec, v_init, dt, T, alpha);

                for (int i = 0; i < n; i++)
                {
                    vec[i] /= H.w_Degree(i);
                }

                int[] index = Enumerable.Range(0, n).ToArray<int>();
                Array.Sort<int>(index, (a, b) => vec[a].CompareTo(vec[b]));

                Array.Reverse(index);

                double vol_V = 0;
                for (int i = 0; i < n; i++) vol_V += H.w_Degree(i);

                var num_contained_nodes = new Dictionary<int, int>();
                for (int eid = 0; eid < H.m; eid++)
                {
                    num_contained_nodes.Add(eid, 0);
                }

                double cut_val = 0;
                double vol_S = 0;
                double conductance = double.MaxValue;
                int best_index = -1;

                foreach (int i in index)
                {
                    vol_S += H.w_Degree(i);
                    if (vol_S <= vol_V / 10.0)
                    {
                        foreach (var e in H.incident_edges[i])
                        {
                            if (num_contained_nodes[e] == 0)
                            {
                                cut_val += H.weights[e];
                            }
                            if (num_contained_nodes[e] == edge_size[e] - 1)
                            {
                                cut_val -= H.weights[e];
                            }
                            num_contained_nodes[e] += 1;
                        }
                        conductance = cut_val / Math.Min(vol_S, vol_V - vol_S);
                        if (conductance < min_conductance)
                        {
                            min_conductance = conductance;
                            best_index = i;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
            time.Stop();
            TimeSpan ts = time.Elapsed;

            Console.WriteLine("conductance: " + min_conductance);
            Console.WriteLine("time(s): " + time.ElapsedMilliseconds/1000.0);
        }



        // for Proposed_local_round
        public static void Proposed_local_round(string fn, int v_init)
        {
            var H = Hypergraph.Open(fn);

            var time = new System.Diagnostics.Stopwatch();
            time.Start();

            int n = H.n;
            int m = H.m;

            const double eps = 0.9;
            
            const double dt = 1.0;
            const double T = 30.0;

            var A_cand = new List<double>();
            for (int i = 0; i <= Math.Log(n * m) / Math.Log(1 + eps); i++)
            {
                A_cand.Add(Math.Pow(1 + eps, i) / (n * m));
            }

            var edge_size = new Dictionary<int, int>();
            for (int eid = 0; eid < H.m; eid++)
            {
                edge_size.Add(eid, H.ID_rev[eid].Count());
            }

            double min_conductance = double.MaxValue;

            foreach (double alpha in A_cand)
            {

                var vec = CreateVector.Dense<double>(n);

                vec[v_init] = 1.0;

                vec = Hypergraph.Simulate_round(H, vec, v_init, dt, T, alpha);

                for (int i = 0; i < n; i++)
                {
                    vec[i] /= H.w_Degree(i);
                }

                int[] index = Enumerable.Range(0, n).ToArray<int>();
                Array.Sort<int>(index, (a, b) => vec[a].CompareTo(vec[b]));

                Array.Reverse(index);

                double vol_V = 0;
                for (int i = 0; i < n; i++) vol_V += H.w_Degree(i);

                var num_contained_nodes = new Dictionary<int, int>();
                for (int eid = 0; eid < H.m; eid++)
                {
                    num_contained_nodes.Add(eid, 0);
                }

                double cut_val = 0;
                double vol_S = 0;
                double conductance = double.MaxValue;
                int best_index = -1;

                foreach (int i in index)
                {
                    vol_S += H.w_Degree(i);
                    if (vol_S <= vol_V / 10.0)
                    {
                        foreach (var e in H.incident_edges[i])
                        {
                            if (num_contained_nodes[e] == 0)
                            {
                                cut_val += H.weights[e];
                            }
                            if (num_contained_nodes[e] == edge_size[e] - 1)
                            {
                                cut_val -= H.weights[e];
                            }
                            num_contained_nodes[e] += 1;
                        }
                        conductance = cut_val / Math.Min(vol_S, vol_V - vol_S);
                        //Console.WriteLine($"{cut_val}, {vol_S}, {vol_V}, {conductance}");
                        if (conductance < min_conductance)
                        {
                            min_conductance = conductance;
                            best_index = i;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
            time.Stop();
            TimeSpan ts = time.Elapsed;

            Console.WriteLine("conductance: " + min_conductance);
            Console.WriteLine("time(s): " + time.ElapsedMilliseconds/1000.0);
        }


        public static void Main(string[] args)
        {
            //Proposed_local("../../instance/dbpedia-writer_LCC.txt", 0);
            Proposed_local_round("../../instance/dbpedia-writer_LCC.txt", 0); // Run Algorithm 1 for DBpedia Writers with seed node 0

        }


    }
}
