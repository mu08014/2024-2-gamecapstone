using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;

//Custom Math
public class CMath
{
    public static double PI_4 = 0.785398163397448309616;
    public static double M_PI = 3.14159265358979323846;
    public static double M_PI_4 = 0.785398163397448309616;
    public static double epsilon = 1e-16;

    public static double cube(double x)
    {
        return x * x * x;
    }

    public static double lerp(double v0, double v1, double t)
    {
        return v0 + t * (v1 - v0);
    }

    public static double clamp(double a, double lower, double upper)
    {
        if (a < lower) return lower;
        else if (a > upper) return upper;
        else return a;
    }

    public static int bipart_closest(in Vectors v, double val)
    {
        List<double> buff = new List<double>(v.size());
        for (int i = 0; i < v.size(); i++)
        {
            buff[i] = v[i];
        }
        buff.Sort();

        int low = 0;
        int high = buff.Count;

        while (low < high)
        {
            int mid = (low + high) / 2;

            if (buff[mid] < val)
                low = mid + 1;
            else
                high = mid;
        }

        if (low >= v.size())
        {
            return v.size() - 1;
        }
        return low;
    }

    public static MatrixXs matProduct(Vectors a, Vectors b)
    {
        int aSize = a.Size;
        int bSize = b.Size;
        MatrixXs result = new MatrixXs(aSize, bSize);
        for (int i = 0; i < aSize; i++)
        {
            for (int j = 0; j < bSize; j++)
            {
                result[i, j] = a[i] * b[j];
            }
        }
        return result;
    }

    public static MatrixXs identity(int size)
    {
        MatrixXs result = new MatrixXs(size, size);
        for (int i = 0; i < size; i++)
        {
            result[i, i] = 1;
        }
        return result;
    }

    public static MatrixXs outerProd(Vectors a, Vectors b)
    {
        MatrixXs result = new MatrixXs(a.Size, b.Size);
        for (int i = 0; i < a.Size; i++)
        {
            for (int j = 0; j < b.Size; j++)
            {
                result[i, j] = a[i] * b[j];
            }
        }
        return result;
    }

    public static MatrixXs crossMat(Vectors a)
    {
        MatrixXs result = new MatrixXs(3, 3);
        result[0, 0] = 0;
        result[0, 1] = -a[2];
        result[0, 2] = a[1];
        result[1, 0] = a[2];
        result[1, 1] = 0;
        result[1, 2] = -a[0];
        result[2, 0] = -a[1];
        result[2, 1] = a[0];
        result[2, 2] = 0;
        return result;
    }

    public static bool isSmall(double x)
    {
        return x < 1e-12;
    }

    public static Vectors orthonormalParallelTransport(in Vectors u, in Vectors t0,
                                  in Vectors t1)
    {
        // This should be called only to transport an orthogonal vector

        Vectors b = t0.cross(t1);
        double bNorm = b.norm();
        if (isSmall(bNorm)) return u;
        b /= bNorm;

        Vectors n0 = t0.cross(b);
        Vectors n1 = t1.cross(b);

        return u.dot(n0) * n1 + u.dot(b) * b;
    }

    public static void rotateAxisAngle(ref Vectors v, in Vectors z, double theta)
    {
        if (theta == 0) return;

        double c = Math.Cos(theta);
        double s = Math.Sin(theta);

        v = c * v + s * z.cross(z) + z.dot(z) * (1.0 - c) * z;
    }

    public static double signedAngle(in Vectors u, in Vectors v, in Vectors n)
    {
        Vectors w = u.cross(v);
        double angle = Math.Atan2(w.norm(), u.dot(v));
        if (n.dot(w) < 0) return -angle;
        return angle;
    }

    public static void orthoNormalize(ref Vectors v, in Vectors n) 
    {
        Vectors vs = v;
        vs -= vs.dot(n)* n;
        vs = vs.normalize();
        v = vs;
    }

    public static Vectors findNormal(Vectors v)
    {
        Vectors result = new Vectors(v.Size);

        int maxCoordinate = 0;
        int n = v.Size;
        for (int i = 0; i < n; i++)
        {
            if (Math.Abs(v[i]) > Math.Abs(v[maxCoordinate]))
            {
                maxCoordinate = i;
            }
        }

        {
            int otherCoordinate = (maxCoordinate + 1) % n;
            result[otherCoordinate] = v[maxCoordinate];
            result[maxCoordinate] = -v[otherCoordinate];
        }

        return result.normalized();
    }

    public static double innerBProduct(in MatrixXs B, in Vectors u, in  Vectors v)
    {
        return B[0, 0] * u[0] * v[0] + B[0, 1] * (u[0] * v[1] + u[1] * v[0]) +
            B[1, 1] * u[1] * v[1];
    }

    public static void symBProduct(int n, MatrixXs result, in MatrixXs B,
                        in MatrixXs Q)
    {
        for (int i = 0; i < n; ++i) {
            Vectors Qrow_i = Q.row(i);
            result[i, i] = innerBProduct(B, Qrow_i, Qrow_i);
            for (int j = 0; j < i; ++j)
                result[i, j] = result[j, i] = innerBProduct(B, Qrow_i, Q.row(j));
        }
    }

    public static void compute_cwiseProduct(ref Vectors b, in Vectors x, in Vectors y)
    {
        Vectors b_ptr = b;
        Vectors x_ptr = x;
        Vectors y_ptr = y;

        Parallel.For(
            0, x.size(), k => { b_ptr[k] = x_ptr[k] * y_ptr[k]; });

        b = b_ptr;
    }

    public static void accumulateJTPhi_coeff(ref Vectors b, in double c, in Vectors x, in SparseXs JT)
    {
        Vectors b_ptr = b;
        Vectors x_ptr = x;
        double c_ptr = c;
        SparseXs JT_ptr = JT;

        Parallel.For(0, JT.outerSize(), k => {
            double sum = 0.0;
            foreach (var ((row, col), val) in JT_ptr.matrix)
            {
                sum += x_ptr[row] * val;
            }
            b_ptr[k] += sum * c_ptr;
        });
    }

    public static void computeJTPhi_coeff(ref Vectors b, in Vectors c,
        in Vectors x, in SparseXs JT)
    {
        Vectors b_ptr = b;
        Vectors c_ptr = c;
        Vectors x_ptr = x;
        SparseXs JT_ptr = JT;

        Parallel.For(0, JT.outerSize(), k => {
            double sum = 0.0;
            foreach (var ((row, col), val) in JT_ptr.matrix)
            {
                sum += x_ptr[row] * val;
            }
            b_ptr[k] = sum * c_ptr[k];
        });
    }
}

[Serializable]
public class Vectors
{
    [SerializeField]
    public double[] elements = new double[3];

    public int Size => elements.Length;

    public static explicit operator UnityEngine.Vector3(Vectors v)
    {
        return new UnityEngine.Vector3((float)v[0], (float)v[1], (float)v[2]);
    }
    public static explicit operator Vectors(UnityEngine.Vector3 v)
    {

        Vectors V = new Vectors(3);
        V.elements[0] = v.x;
        V.elements[1] = v.y;
        V.elements[2] = v.z;
        return V;
    }

    public Vectors()
    {
        elements = new double[3];
    }

    public Vectors(int size)
    {
        elements = new double[size];
    }

    public Vectors(double[] list)
    {
        elements = list;
    }

    public Vectors(double[] list, int start, int size)
    {
        if (list.Length == size)
            elements = list;
        else
        {
            elements = new double[list.Length - size];
            for (int i = 0; i < list.Length - size; i++)
            {
                elements[i] = (list[start + i]);
            }
        }
    }

    public Vectors(List<double> list, int start, int size)
    {
        if (list.Count == size)
            elements = list.ToArray();
        else
        {
            elements = new double[list.Count - size];
            for (int i = 0; i < list.Count - size; i++)
            {
                elements[i] = (list[start + i]);
            }
        }
    }

    public double this[int index]
    {
        get => elements[index];
        set => elements[index] = value;
    }

    public static Vectors operator +(Vectors a, Vectors b)
    {
        if (a.Size != b.Size)
            throw new System.Exception("벡터의 크기가 다릅니다.");

        Vectors result = new Vectors(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    public static Vectors operator -(Vectors a, Vectors b)
    {
        if (a.Size != b.Size)
            throw new System.Exception("벡터의 크기가 다릅니다.");

        Vectors result = new Vectors(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    public static Vectors operator *(double a, Vectors b)
    {
        Vectors result = new Vectors(b.Size);
        for (int i = 0; i < b.Size; i++)
        {
            result[i] = a * b[i];
        }
        return result;
    }

    public static Vectors operator *(Vectors b, double a)
    {
        Vectors result = new Vectors(b.Size);
        for (int i = 0; i < b.Size; i++)
        {
            result[i] = a * b[i];
        }
        return result;
    }

    public static Vectors operator *(Vectors a, Vectors b)
    {
        if (a.Size != 4 || b.Size != 4)
            throw new System.Exception("쿼터니언 곱을 할 수 없습니다.");
        Vectors result = new Vectors(a.Size);
        result[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
        result[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];
        result[2] = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1];
        result[3] = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0];

        return result;
    }

    public static Vectors operator /(Vectors a, double b)
    {
        Vectors result = new Vectors(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            if (a[i] != 0)
                result[i] = a[i] / b;
        }
        return result;
    }

    public static Vectors operator +(UnityEngine.Quaternion q, Vectors v)
    {
        Vectors result = new Vectors(3);
        for (int i = 0;i < 3;i++)
        {
            result[0] = q.x + v[0];
            result[1] = q.y + v[1];
            result[2] = q.z + v[2];
        }
        return result;
    }

    public void setZero()
    {
        elements = new double[elements.Length];
    }

    public int size()
    {
        return elements.Length;
    }

    public double norm()
    {
        double result = 0;
        for (int i = 0; i < Size; i++)
        {
            result += elements[i] * elements[i];
        }
        result = Math.Sqrt(result);
        return result;
    }

    public Vectors qinverse()
    {
        if (Size != 4)
            throw new System.Exception("쿼터니언이 아닙니다.");

        Vectors result = new Vectors(4);

        double qsquare = norm() * norm();

        if (qsquare != 0)
        {
            for (int i = 0; i < 4; i++)
            {
                if (elements[i] != 0)
                {
                    result[i] = elements[0] / qsquare;
                }

                if (i != 0)
                {
                    result[i] *= -1;
                }
            }   
        }

        return result;
    }

    public void resize(int size)
    {
        elements = new double[size];
    }

    public Vectors segment(int start, int n)
    {
        Vectors vectors = new Vectors(n);
        for (int i = start; i < start + n; i++)
        {
            vectors[i - start] = elements[i];
        }
        return vectors;
    }

    public void SetSegment(int start, int n, Vectors v)
    {
        if (v.Size != n)
            throw new System.Exception("벡터의 크기가 다릅니다.");

        for (int i = start; i < start + n; i++)
        {
            if (i >= elements.Length)
            {
                throw new System.Exception("벡터의 크기를 초과했습니다.");
            }
            else
            {
                elements[i] = v[i - start];
            }
        }
    }

    public Vectors normalize()
    {
        double sum = 0;
        Vectors vectors = new Vectors(Size);
        for (int i = 0; i < Size; i++)
        {
            sum += elements[i] * elements[i];
        }

        if (sum > 0)
        {
            sum = Math.Sqrt(sum);
            for (int i = 0; i < Size; i++)
            {
                if (elements[i] != 0)
                    vectors[i] = elements[i] / sum;
            }
        }
        return vectors;
    }

    public Vectors normalized()
    {
        double sum = 0;
        Vectors vectors = new Vectors(Size);
        for (int i = 0; i < Size; i++)
        {
            sum += elements[i] * elements[i];
        }

        if (sum > 0)
        {
            sum = Math.Sqrt(sum);
            for (int i = 0; i < Size; i++)
            {
                if (elements[i] != 0)
                    vectors[i] = elements[i] / sum;
            }
        }
        return vectors;
    }

    public Vectors cross(Vectors v)
    {
        if (v.Size != 3)
            throw new System.Exception("벡터의 크기가 다릅니다.");
        Vectors result = new Vectors(3);
        result[0] = elements[1] * v[2] - elements[2] * v[1];
        result[1] = elements[2] * v[0] - elements[0] * v[2];
        result[2] = elements[0] * v[1] - elements[1] * v[0];
        return result;
    }

    public Vectors setConstant(double value)
    {
        for (int i = 0; i < Size; i++)
        {
            elements[i] = value;
        }
        return this;
    }

    public double squaredNorm()
    {
        double sum = 0;
        for (int i = 0; i < Size; i++)
        {
            sum += elements[i] * elements[i];
        }
        return sum;
    }

    public Vectors Clone()
    {
        Vectors result = new Vectors(Size);
        for (int i = 0; i < Size; i++)
        {
            result[i] = elements[i];
        }

        return result;
    }

    public double dot(Vectors a)
    {
        if (a.Size != elements.Length)
            throw new System.Exception("벡터의 크기가 다릅니다.");

        double result = 0;
        for (int i = 0; i < Size; i++)
        {
            result += elements[i] * a[i];
        }
        return result;
    }

    public Vectors cwiseInverse()
    {
        Vectors result = new Vectors(Size);
        for (int i = 0; i < Size;i++)
        {
            result[i] = elements[i] != 0 ? 1.0 / elements[i] : 0;
        }
        return result;
    }

    public Vectors cwiseProduct(Vectors v)
    {
        if (v.Size != elements.Length)
            throw new System.Exception("벡터의 크기가 다릅니다.");

        Vectors result = new Vectors(Size);
        for(int i = 0; i < Size; i++)
        {
            result[i] = v[i] * elements[i];
        }
        return result;
    }

    public void PlusSegment(int start, int n, Vectors v)
    {
        if (v.Size != n)
            throw new System.Exception("벡터의 크기가 다릅니다.");

        for (int i = start; i < start + n; i++)
        {
            if (i >= elements.Length)
            {
                throw new System.Exception("벡터의 크기를 초과했습니다.");
            }
            else
            {
                elements[i] += v[i - start];
            }
        }
    }

    public void conservativeResize(int size)
    {
        int currentSize = elements.Length;

        if (size > currentSize)
        {
            // 크기가 증가하는 경우, 부족한 부분을 0으로 채움
            double[] newArray = new double[size];
            System.Array.Copy(elements, 0, newArray, 0, currentSize);

            elements = newArray;
        }
        else if (size < currentSize)
        {
            // 크기가 감소하는 경우, 요소를 잘라냄
            double[] newArray = new double[size];

            System.Array.Copy(elements, 0, newArray, 0, size);
            elements = newArray;
        }
    }
}

public class VectorXi
{
    public List<int> elements;
    public int Size => elements.Count;

    public VectorXi(int size)
    {
        elements = new List<int>(new int[size]);
    }

    public VectorXi(int[] list)
    {
        elements = new List<int>(list);
    }

    public int this[int index]
    {
        get => elements[index];
        set => elements[index] = value;
    }

    public static VectorXi operator +(VectorXi a, VectorXi b)
    {
        if (a.Size != b.Size)
            throw new System.Exception("벡터의 크기가 다릅니다.");

        VectorXi result = new VectorXi(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    public static VectorXi operator -(VectorXi a, VectorXi b)
    {
        if (a.Size != b.Size)
            throw new System.Exception("벡터의 크기가 다릅니다.");

        VectorXi result = new VectorXi(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    public void SetSegment(int start, int n, VectorXi v)
    {
        if (v.Size != n)
            throw new System.Exception("벡터의 크기가 다릅니다.");
        for (int i = start; i < start + n; i++)
        {
            if (i >= elements.Count)
            {
                elements.Add(v[i - start]);
            }
            else
            {
                elements[i] = v[i - start];
            }
        }
    }

    public VectorXi segment(int start, int n)
    {
        VectorXi vectors = new VectorXi(n);
        for (int i = start; i < start + n; i++)
        {
            vectors[i - start] = elements[i];
            
        }
        return vectors;
    }

    public VectorXi setConstant(int value)
    {
        for (int i = 0; i < Size; i++)
        {
            elements[i] = value;
        }
        return this;
    }

    public void resize(int size)
    {
        elements = new List<int>(new int[size]);
    }

    public void conservativeResize(int size)
    {
        int currentSize = elements.Count;

        if (size > currentSize)
        {
            // 크기가 증가하는 경우, 부족한 부분을 0으로 채움
            elements.AddRange(new int[size - currentSize]);
        }
        else if (size < currentSize)
        {
            // 크기가 감소하는 경우, 요소를 잘라냄
            elements.RemoveRange(size, currentSize - size);
        }
    }

    public VectorXi Clone()
    {
        VectorXi result = new VectorXi(Size);
        for (int i = 0; i < Size; i++)
        {
            result[i] = elements[i];
        }

        return result;
    }
}

//희소 행렬에 사용
public class Triplets
{
    private int row { get; set; }
    private int col { get; set; }
    private double value { get; set; }

    public Triplets(int r, int c, double v)
    {
        if (r < 0 || c < 0)
            Debug.LogError("Triplets 의 생성하려는 row나 col 이 0보다 작음");
        row = r;
        col = c;
        value = v;
    }

    public int getrow()
    {
        return row;
    }

    public int getcol()
    {
        return col;
    }

    public double getvalue()
    {
        return value;
    }
}

public class TripletXs
{
    public List<Triplets> elements;

    public int Size => elements.Count;

    public TripletXs()
    {
        elements = new List<Triplets>();
    }

    public TripletXs(int size)
    {
        elements = new List<Triplets>(size);
        for (int i = 0; i < size; i++)
        {
            elements.Add(new Triplets(0, 0, 0));
        }
    }

    public Triplets this[int index]
    {
        get => elements[index];
        set => elements[index] = value;
    }

    public void Add(int row, int col, int value)
    {
        Triplets triplets = new Triplets(row, col, value);
        elements.Add(triplets);
    }

    public void Add(Triplets triplets)
    {
        elements.Add(triplets);
    }

    public void resize(int size)
    {
        int currentSize = elements.Count;

        if (size > currentSize)
        {
            elements.AddRange(new int[size - currentSize]);
        }
        else if (size < currentSize)
        {
            elements.RemoveRange(size, currentSize - size);
        }
    }

    public void reserve(int size)
    {
        elements.Capacity = size;
    }
}

public class MatrixXs
{
    private int rows; // 행 개수
    private int cols; // 열 개수
    private double[,] data; // 행렬 데이터

    public MatrixXs()
    {
        rows = 0;
        cols = 0;
        data = new double[rows, cols];
    }

    public MatrixXs(int rows, int cols)
    {
        this.rows = rows;
        this.cols = cols;
        data = new double[rows, cols];
    }

    public MatrixXs(MatrixXs other)
    {
        rows = other.rows;
        cols = other.cols;
        data = new double[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                data[i, j] = other.data[i, j];
    }

    public static MatrixXs operator*(double b, MatrixXs a)
    {
        MatrixXs result = new MatrixXs(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
            {
                result[i, j] = a.data[i, j] * b;
            }
        }
        return result;
    }

    public static MatrixXs operator *(MatrixXs a, double b)
    {
        MatrixXs result = new MatrixXs(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
            {
                result[i, j] = a.data[i, j] * b;
            }
        }
        return result;
    }

    public static MatrixXs operator /(MatrixXs a, double b)
    {
        if (b == 0)
            Debug.LogError("행렬을 0으로 나눌 수 없습니다");

        MatrixXs result = new MatrixXs(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
            {
                if (a.data[i, j] != 0)
                    result[i, j] = a.data[i, j] / b;
            }
        }
        return result;
    }

    public static MatrixXs operator +(MatrixXs a, MatrixXs b)
    {
        if (a.rows != b.rows || a.cols != b.cols)
        {
            Debug.LogError("행렬 덧셈을 하기 위해서 두 행렬의 크기가 같아야 합니다.");
        }
        MatrixXs result = new MatrixXs(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
            {
                result[i, j] = a.data[i, j] + b.data[i, j];
            }
        }
        return result;
    }

    public static MatrixXs operator -(MatrixXs a, MatrixXs b)
    {
        if (a.rows != b.rows || a.cols != b.cols)
        {
            Debug.LogError("행렬 셈을 하기 위해서 두 행렬의 크기가 같아야 합니다.");
        }
        MatrixXs result = new MatrixXs(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
            {
                result[i, j] = a.data[i, j] - b.data[i, j];
            }
        }
        return result;
    }

    public void Setblock(int row_size, int col_size, int start_row, int start_col, MatrixXs m)
    {
        for (int i = start_row; i < start_row + row_size; i++)
        {
            for (int j = start_col; j < start_col + col_size; j++)
            {
                data[i, j] = m[i - start_row, j - start_col];
            }
        }
    }

    public void Setblock(int row_size, int col_size, int start_row, int start_col, Vectors v)
    {
        if (col_size != 1)
        {
            Debug.LogError("열 크기가 1이 아니면, 벡터로 값 설정할 수 없습니다.");
        }

        for (int i = 0; i < row_size; i++)
        {
            data[i + start_row, start_col] = v[i];
        }
    }

    public MatrixXs block(int row_size, int col_size, int start_row, int start_col)
    {
        MatrixXs result = new MatrixXs(row_size, col_size);
        for (int i = start_row; i < start_row + row_size; i++)
        {
            for (int j = start_col; j < start_col + col_size; j++)
            {
                result[i - start_row, j - start_col] = data[i, j];
            }
        }
        return result;
    }

    public MatrixXs Clone()
    {
        return new MatrixXs(this);
    }

    public void SetElement(int row, int col, double value)
    {
        if (row >= 0 && row < rows && col >= 0 && col < cols)
        {
            data[row, col] = value;
        }
        else
        {
            throw new IndexOutOfRangeException("인덱스가 범위를 벗어났습니다.");
        }
    }

    // 요소를 가져오는 메서드
    public double GetElement(int row, int col)
    {
        if (row >= 0 && row < rows && col >= 0 && col < cols)
        {
            return data[row, col];
        }
        else
        {
            throw new IndexOutOfRangeException("인덱스가 범위를 벗어났습니다.");
        }
    }

    public int GetRows()
    {
        return rows;
    }

    public int GetCols()
    {
        return cols;
    }

    public void setZero()
    {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                data[i, j] = 0;
    }


    // 행렬 곱셈 메서드
    public static MatrixXs operator *(MatrixXs a, MatrixXs b)
    {
        return MatrixMultiplicationGPU.Multiply(a, b);
    }

    public double this[int i, int j]
    {
        get => data[i, j];
        set => data[i, j] = value;
    }

    public MatrixXs transpose()
    {
        MatrixXs matrix = new MatrixXs(cols, rows);
        for (int i = 0; i < cols; i++)
            for (int j = 0; j < rows; j++)
                matrix.data[i, j] = data[j, i];
        return matrix;
    }

    public Vectors row(int i)
    {
        Vectors v = new Vectors(cols);
        for (int j = 0; j < cols; j++)
        {
            v[j] = data[i, j];
        }
        return v;
    }

    public Vectors col(int i)
    {
        Vectors v = new Vectors(i);
        for (int j = 0; j < rows; j++)
        {
            v[j] = data[j, i];
        }
        return v;
    }

    public void resize(int rows, int cols)
    {
        this.rows = rows;
        this.cols = cols;
        data = new double[rows, cols];
        setZero();
    }
}

public class Rotation2D
{
    private double angle;  // 회전 각도 (라디안)

    public Rotation2D(double angleInRadians)
    {
        angle = angleInRadians;
    }

    // 회전 행렬을 사용하여 벡터 (x, y)를 회전
    public (double x, double y) Rotate(double x, double y)
    {
        double cosTheta = Math.Cos(angle);
        double sinTheta = Math.Sin(angle);

        double newX = cosTheta * x - sinTheta * y;
        double newY = sinTheta * x + cosTheta * y;

        return (newX, newY);
    }

    public MatrixXs toRotationMatrix()
    {
        double cosTheta = Math.Cos(angle);
        double sinTheta = Math.Sin(angle);
        MatrixXs matrix = new MatrixXs(2, 2);

        matrix[0, 0] = cosTheta;
        matrix[0, 1] = -sinTheta;
        matrix[1, 0] = sinTheta;
        matrix[1, 1] = cosTheta;

        return matrix;
    }
}

public class SparseXs
{
    public Dictionary<(int, int), double> matrix;
    public int Rows { get; private set; }
    public int Cols { get; private set; }

    public SparseXs(int rows, int cols)
    {
        matrix = new Dictionary<(int, int), double>();
        Rows = rows;
        Cols = cols;
    }

    public SparseXs(SparseXs spXs)
    {
        matrix = new Dictionary<(int, int), double>();
        foreach (var d in spXs.matrix)
        {
            if (Math.Abs(d.Value) > CMath.epsilon)
                matrix.Add((d.Key.Item1, d.Key.Item2), d.Value);
        }
        Rows = spXs.Rows;
        Cols = spXs.Cols;
    }

    public static SparseXs operator +(SparseXs a, SparseXs b)
    {
        if (a.Rows != b.Rows || a.Cols != b.Cols)
            throw new ArgumentException("행렬 덧셈을 위해선 A의 행열 개수가 B의 행열 개수와 같아야 합니다.");

        SparseXs result = new SparseXs(a.Rows, b.Cols);

        foreach (var ((rowA, colA), valueA) in a.matrix)
        {
            result.insert(rowA, colA, valueA);
        }

        foreach (var ((rowB, colB), valueB) in b.matrix)
        {
            if (result.matrix.ContainsKey((rowB, colB)))
            {
                if (Math.Abs(result.matrix[(rowB, colB)] + valueB) <= CMath.epsilon)
                    result.matrix.Remove((rowB, colB));
                else
                    result.matrix[(rowB, colB)] += valueB;
            }
            else
            {
                result.insert(rowB, colB, valueB);
            }
        }

        return result;
    }

    public static SparseXs operator *(SparseXs b, double a)
    {
        SparseXs result = new SparseXs(b.Rows, b.Cols);
        if (Math.Abs(a) > CMath.epsilon)
        {
            foreach (var ((row, col), value) in b.matrix)
            {
                result.insert(row, col, value * a);
            }
        }
        return result;
    }

    public static SparseXs operator *(double a, SparseXs b)
    {
        SparseXs result = new SparseXs(b.Rows, b.Cols);
        if (Math.Abs(a) > CMath.epsilon)
        {
            foreach (var ((row, col), value) in b.matrix)
            {
                result.insert(row, col, value * a);
            }
        }
        return result;
    }

    public static Vectors operator *(SparseXs a, Vectors b)
    {
        if (a.Cols != b.Size)
            throw new ArgumentException("행렬 곱셈을 위해선 A의 열 개수가 B의 벡터 크기가 같아야 합니다.");

        Vectors result = new Vectors(a.Rows);
        foreach (var ((rowA, colA), valueA) in a.matrix)
        {
            if (Math.Abs(valueA) > CMath.epsilon)
                result[rowA] += b[colA] * valueA;
        }
        return result;
    }

    public static SparseXs operator *(SparseXs a, SparseXs b)
    {
        if (a.Cols != b.Rows)
            throw new ArgumentException("행렬 곱셈을 위해선 A의 열 개수가 B의 행 개수와 같아야 합니다.");

        SparseXs result = new SparseXs(a.Rows, b.Cols);

        var bRowDict = new Dictionary<int, List<(int col, double value)>>();

        foreach (var ((rowB, colB), valueB) in b.matrix)
        {
            if (!bRowDict.ContainsKey(rowB))
                bRowDict[rowB] = new List<(int col, double value)>();

            bRowDict[rowB].Add((colB, valueB));
        }

        foreach (var ((rowA, colA), valueA) in a.matrix)
        {
            if (!bRowDict.ContainsKey(colA)) continue;

            foreach (var (colB, valueB) in bRowDict[colA])
            {
                // rowA와 colB에 대해 결과 계산
                double newValue = result.GetValue(rowA, colB) + valueA * valueB;
                result.insert(rowA, colB, newValue);
            }
        }

        return result;
    }

    public void insert(int row, int col, double value)
    {
        if (row >= Rows || col >= Cols || row < 0 || col < 0)
        {
            Debug.LogError($"Invalid index: row={row}, col={col}");
        }
        if (Math.Abs(value) > CMath.epsilon)
            matrix[(row, col)] = value;
    }

    public double GetValue(int row, int col)
    {
        if (row >= Rows || col >= Cols || row < 0 || col < 0)
        {
            Debug.LogError($"Invalid index: row={row}, col={col}");
        }

        return matrix.ContainsKey((row, col)) ? matrix[(row, col)] : 0.0;
    }

    public void setFromTriplets(in TripletXs tpXs)
    {
        for (int i = 0; i < tpXs.Size; i++)
        {
            if (Math.Abs(tpXs[i].getvalue()) > CMath.epsilon)
            {
                var key = (tpXs[i].getrow(), tpXs[i].getcol());
                matrix[key] = tpXs[i].getvalue();
            }
        }

    }

    public void resize(int row, int col)
    {
        Rows = row; Cols = col;
        matrix = new Dictionary<(int, int), double>();
    }

    public SparseXs transpose()
    {
        SparseXs tp = new SparseXs(Cols, Rows);
        foreach (var ((row, col), value) in matrix)
        {
            tp.insert(col, row, value);
        }

        return tp;
    }

    public Vectors diagonal()
    {
        Vectors result = new Vectors(Rows > Cols ? Cols : Rows);
        for (int i = 0; i < Math.Min(Rows, Cols); i++)
        {
            result[i] = GetValue(i, i);
        }

        return result;
    }

    public int outerSize()
    {
        return Cols < Rows ? Cols : Rows;
    }

}

public static class MatrixMultiplicationGPU
{
    public static ComputeShader computeShader;

    // 행렬 곱셈을 위한 GPU 연산 실행
    public static MatrixXs Multiply(MatrixXs a, MatrixXs b)
    {
        computeShader = Resources.Load<ComputeShader>("MatrixMultiplication");

        if (a.GetCols() != b.GetRows())
        {
            throw new InvalidOperationException("첫 번째 행렬의 열의 개수와 두 번째 행렬의 행의 개수가 일치해야 합니다.");
        }

        int rows = a.GetRows();
        int cols = b.GetCols();
        int common = a.GetCols();

        MatrixXs result = new MatrixXs(rows, cols);

        float[] aData = new float[rows * common];
        float[] bData = new float[common * cols];
        float[] resultData = new float[rows * cols];

        for (int i = 0; i < rows; i++)
        {
            for (int k = 0; k < common; k++)
            {
                aData[i * common + k] = (float)a.GetElement(i, k);
            }
        }

        for (int k = 0; k < common; k++)
        {
            for (int j = 0; j < cols; j++)
            {
                bData[k * cols + j] = (float)b.GetElement(k, j);
            }
        }

        ComputeBuffer bufferA = new ComputeBuffer(rows * common, sizeof(float));
        ComputeBuffer bufferB = new ComputeBuffer(common * cols, sizeof(float));
        ComputeBuffer bufferResult = new ComputeBuffer(rows * cols, sizeof(float));

        bufferA.SetData(aData);
        bufferB.SetData(bData);

        int kernelHandle = computeShader.FindKernel("CSMatrixMultiply");
        computeShader.SetBuffer(kernelHandle, "MatrixA", bufferA);
        computeShader.SetBuffer(kernelHandle, "MatrixB", bufferB);
        computeShader.SetBuffer(kernelHandle, "Result", bufferResult);
        computeShader.SetInt("Rows", rows);
        computeShader.SetInt("Cols", cols);
        computeShader.SetInt("Common", common);

        computeShader.Dispatch(kernelHandle, Mathf.CeilToInt((float)cols / 8), Mathf.CeilToInt((float)rows / 8), 1);

        bufferResult.GetData(resultData);

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                result.SetElement(i, j, resultData[i * cols + j]);
            }
        }

        bufferA.Release();
        bufferB.Release();
        bufferResult.Release();

        return result;
    }
}