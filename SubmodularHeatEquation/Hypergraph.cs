using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace SubmodularHeatEquation
{
    public class Hypergraph
    {
        public List<List<int>> edges = new List<List<int>>();
        public List<List<int>> incident_edges = new List<List<int>>();
        public List<double> weights = new List<double>();
        public Dictionary<List<int>, int> ID = new Dictionary<List<int>, int>();
        public Dictionary<int, List<int>> ID_rev = new Dictionary<int, List<int>>();

        public static Hypergraph RandomHypergraph(int n, int m, int k)
        {
            var H = new Hypergraph();

            for (int i = 0; i < m; i++)
            {
                var edge = new List<int>();
                for (int j = 0; j < k; j++)
                {
                    int v = (int)(Util.Xor128() % n);
                    edge.Add(v);
                }
                edge = edge.Distinct().ToList();
                H.AddEdge(edge);
            }
            return H;
        }

        public int Degree(int v)
        {
            return incident_edges[v].Count;
        }

        public double w_Degree(int v)
        {
            double sum = 0;
            foreach (int e in incident_edges[v])
            {
                sum += weights[e];
            }
            return sum;
        }

        public int n
        {
            get
            {
                return incident_edges.Count;
            }
        }

        public int m
        {
            get
            {
                return edges.Count;
            }
        }

        public void AddEdge(List<int> edge, double w = 1)
        {
            int eid = edges.Count;
            edges.Add(edge);
            foreach (var v in edge)
            {
                while (v >= incident_edges.Count)
                {
                    incident_edges.Add(new List<int>());
                }
                incident_edges[v].Add(eid);
            }
            weights[eid] = w;
            ID[edge] = eid;
        }

        public static Hypergraph Open(string fn)
        {
            var fs = new FileStream(fn, FileMode.Open);
            var sr = new StreamReader(fs);

            var edges = new List<List<int>>();
            var weights = new List<double>();
            int vertex_num = 0;
            for (string line; (line = sr.ReadLine()) != null;)
            {
                var words = line.Split();
                var edge = new List<int>();
                int i = 0;
                foreach (var word in words)
                {
                    if (i < words.Length - 1)
                    {
                        int v = int.Parse(word);
                        edge.Add(v);
                        if (v >= vertex_num) vertex_num = v + 1;
                        i++;
                    }
                }
                edges.Add(edge);
                weights.Add(double.Parse(words.Last()));
            }

            var H = new Hypergraph();
            H.edges = edges;
            H.weights = weights;
            for (int v = 0; v < vertex_num; v++)
            {
                H.incident_edges.Add(new List<int>());
            }
            for (int i = 0; i < edges.Count; i++)
            {
                var edge = edges[i];
                foreach (var v in edge)
                {
                    H.incident_edges[v].Add(i);
                }
                H.ID.Add(edge, i);
                H.ID_rev.Add(i, edge);
            }
            fs.Close();
            return H;
        }

        public Vector<double> T(Vector<double> vec, int v_init, double alpha)
        {
            var res = CreateVector.Dense<double>(n);

            const double eps = 1e-8;
            foreach (var edge in edges)
            {
                var argmaxs = new List<int>();
                var argmins = new List<int>();
                double maxval = double.MinValue, minval = double.MaxValue;
                foreach (var v in edge)
                {
                    var val = vec[v] / w_Degree(v);
                    if (val > maxval + eps)
                    {
                        maxval = val;
                        argmaxs.Clear();
                        argmaxs.Add(v);
                    }
                    else if (val > maxval - eps)
                    {
                        argmaxs.Add(v);
                    }

                    if (val < minval - eps)
                    {
                        minval = val;
                        argmins.Clear();
                        argmins.Add(v);
                    }
                    else if (val < minval + eps)
                    {
                        argmins.Add(v);
                    }
                }
                foreach (var v in argmaxs)
                {
                    res[v] += weights[ID[edge]] * (maxval - minval) / argmaxs.Count;
                }
                foreach (var v in argmins)
                {
                    res[v] -= weights[ID[edge]] * (maxval - minval) / argmins.Count;
                }
            }

            var res_init = vec;
            res_init[v_init] -= 1;

            var mix = (1 - alpha) * res + alpha * res_init;

            return mix;
        }

        public static Vector<double> Iterate(Hypergraph H, Vector<double> vec, int v_init, double dt, double alpha)
        {
            var dv = H.T(vec, v_init, alpha);
            var res = vec;
            res -= dv * dt;
            return res;
        }

        public static Vector<double> Simulate(Hypergraph H, Vector<double> vec, int v_init, double dt, double T, double alpha)
        {
            var cur_time = 0.0;
            while (cur_time < T)
            {
                var next_time = Math.Min(cur_time + dt, T);
                vec = Iterate(H, vec, v_init, next_time - cur_time, alpha);
                cur_time = next_time;
            }
            return vec;
        }






        public Vector<double> T_round(Vector<double> vec, int v_init, double alpha, List<List<int>> active_edges)
        {
            var res = CreateVector.Dense<double>(n);

            const double eps = 1e-8;

            foreach (var edge in active_edges)
            {
                var argmaxs = new List<int>();
                var argmins = new List<int>();
                double maxval = double.MinValue, minval = double.MaxValue;
                foreach (var v in edge)
                {
                    var val = vec[v] / w_Degree(v);
                    if (val > maxval + eps)
                    {
                        maxval = val;
                        argmaxs.Clear();
                        argmaxs.Add(v);
                    }
                    else if (val > maxval - eps)
                    {
                        argmaxs.Add(v);
                    }

                    if (val < minval - eps)
                    {
                        minval = val;
                        argmins.Clear();
                        argmins.Add(v);
                    }
                    else if (val < minval + eps)
                    {
                        argmins.Add(v);
                    }
                }
                foreach (var v in argmaxs)
                {
                    res[v] += weights[ID[edge]] * (maxval - minval) / argmaxs.Count;
                }
                foreach (var v in argmins)
                {
                    res[v] -= weights[ID[edge]] * (maxval - minval) / argmins.Count;
                }
            }

            var res_init = CreateVector.Dense<double>(n);
            vec.CopyTo(res_init);
            res_init[v_init] -= 1;

            var mix = (1 - alpha) * res + alpha * res_init;

            return mix;
        }

        public static Vector<double> Iterate_round(Hypergraph H, Vector<double> vec, int v_init, double dt, double alpha, List<List<int>> active_edges, List<int> active)
        {
            var dv = H.T_round(vec, v_init, alpha, active_edges);
            var res = vec;
            res -= dv * dt;

            for (int i = 0; i < H.n; i++)
            {
                if (res[i] < 1e-5)
                {
                    res[i] = 0;
                }
            }
            
            var new_active_edges = new List<int>();

            for (int i = 0; i < H.n; i++)
            {
                if (vec[i] == 0 && res[i] != 0)
                {
                    foreach (var f in H.incident_edges[i])
                    {
                        if (active[f] == 0)
                        {
                            new_active_edges.Add(f);
                            active[f] = 1;
                        }
                    }
                }
            }

            foreach (var e in new_active_edges)
            {
                active_edges.Add(H.ID_rev[e]);
            }

            return res;
        }

        public static Vector<double> Simulate_round(Hypergraph H, Vector<double> vec, int v_init, double dt, double T, double alpha)
        {

            var active_edges = new List<List<int>>();
            var active = new List<int>(new int[H.m]);

            foreach (var e in H.incident_edges[v_init])
            {
                active_edges.Add(H.ID_rev[e]);
                active[e] = 1;
            }

            var cur_time = 0.0;
            while (cur_time < T)
            {
                var next_time = Math.Min(cur_time + dt, T);
                vec = Iterate_round(H, vec, v_init, next_time - cur_time, alpha, active_edges, active);
                cur_time = next_time;
            }
            return vec;
        }


        public static double Sqr(double v) { return v * v; }

    }


    public class Graph
    {
        public List<Dictionary<int, double>> adj_list = new List<Dictionary<int, double>>();

        public int Degree(int v)
        {
            return adj_list[v].Count;
        }

        public double w_Degree(int v)
        {
            double sum = 0;
            foreach (var neighbor_val in adj_list[v].Values)
            {
                sum += neighbor_val;
            }
            return sum;
        }

        public int n
        {
            get
            {
                return adj_list.Count;
            }
        }

        public int m
        {
            get
            {
                int sum = 0;
                for (int i = 0; i < n; i++)
                {
                    sum += adj_list[i].Keys.Count;
                }
                return sum / 2;
            }
        }

        public void AddEdge(List<int> edge, double w = 1)
        {
            while (Math.Max(edge[0], edge[1]) >= adj_list.Count)
            {
                adj_list.Add(new Dictionary<int, double>());
            }

            adj_list[edge[0]][edge[1]] = w; 
            adj_list[edge[1]][edge[0]] = w;
        }

    }


}
