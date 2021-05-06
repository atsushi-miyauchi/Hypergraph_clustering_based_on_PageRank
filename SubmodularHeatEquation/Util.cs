using System;
namespace SubmodularHeatEquation
{
    public class Util
    {

        

        public static uint Xor128()
        {
            uint t;
            t = x ^ (x << 11);
            x = y; y = z; z = w;
            return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
        }

        static uint x = 123456789;
        static uint y = 362436069;
        static uint z = 521288629;
        static uint w = 88675123;


        public Util()
        {
        }
    }
}
