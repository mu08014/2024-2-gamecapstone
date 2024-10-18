using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Numerics;
using System.Xml.Linq;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;

//Custom Math
public class CMath
{
    public static double PI_4 = 0.785398163397448309616;
    public static double M_PI = 3.14159265358979323846;
    public static double M_PI_4 = 0.785398163397448309616;

    public static double cube(double x)
    {
        return x * x * x;
    }

    public static double lerp(double v0, double v1, double t)
    {
        return v0 + t * (v1 - v0);
    }

    public static int bipart_closest(in VectorXs v, double val)
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

    public static Vectors orthonormalParallelTransport(in Vectors u, in Vectors t0,
                                  in Vectors t1)
    {
        // This should be called only to transport an orthogonal vector

        Vectors b = t0.cross(t1);
        double bNorm = b.norm();
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

    public static void orthoNormalize(ref VectorXs v, in Vectors n) 
    {
        Vectors vs = v.ToVectors();
        vs -= vs.dot(n)* n;
        vs = vs.normalize();
        v = vs.ToVectorXs();
    }

    public static VectorXs findNormal(VectorXs v)
    {
        VectorXs result = new VectorXs(v.Size);

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

    public static Vectors findNormal(Vectors v)
    {
        Vectors result = new Vectors(v.DIM);
        int maxCoordinate = 0;
        int n = v.DIM;
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
}

public class VectorXs
{
    private List<double> elements;
    public int Size => elements.Count;

    public VectorXs()
    {
        elements = new List<double>();
    }

    public VectorXs(int size)
    {
        elements = new List<double>(new double[size]);
    }

    public double this[int index]
    {
        get => elements[index];
        set => elements[index] = value;
    }

    public static VectorXs operator +(VectorXs a, VectorXs b)
    {
        if (a.Size != b.Size)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");

        VectorXs result = new VectorXs(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    public static VectorXs operator -(VectorXs a, VectorXs b)
    {
        if (a.Size != b.Size)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");

        VectorXs result = new VectorXs(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    public static VectorXs operator *(double a, VectorXs b)
    {
        VectorXs result = new VectorXs(b.Size);
        for (int i = 0; i < b.Size; i++)
        {
            result[i] = a * b[i];
        }
        return result;
    }

    public static VectorXs operator *(VectorXs b, double a)
    {
        VectorXs result = new VectorXs(b.Size);
        for (int i = 0; i < b.Size; i++)
        {
            result[i] = a * b[i];
        }
        return result;
    }

    public static VectorXs operator /(VectorXs a, double b)
    {
        VectorXs result = new VectorXs(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            if (a[i] != 0)
                result[i] = a[i] / b;
        }
        return result;
    }

    public static VectorXs operator +(UnityEngine.Quaternion q, VectorXs v)
    {
        VectorXs result = new VectorXs(3);
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
        for (int i = 0; i < elements.Count; i++)
        {
            elements[i] = 0;
        }
    }

    public int size()
    {
        return elements.Count;
    }

    public void resize(int size)
    {
        List<double> values = new List<double>(Size);
        for (int i = 0; i < Size; i++)
        {
            values.Add(elements[i]);
        }

        elements = new List<double>(size);
        for (int i = 0; i < size; i++)
        {
            if (i < values.Count)
                elements.Add(values[i]);
            else
                elements.Add(0);
        }
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

    // .segment() = Vectors ������ �� ���
    public void SetSegment(int start, int n, Vectors v)
    {
        if (v.DIM != n)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");

        for (int i = start; i < start + n; i++)
        {
            if (i >= elements.Count)
            {
                elements.Add(elements[i]);
            }
            else
            {
                elements[i] = v[i - start];
            }
        }
    }

    public void SetSegment(int start, int n, VectorXs v)
    {
        if (v.Size != n)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");

        for (int i = start; i < start + n; i++)
        {
            if (i >= elements.Count)
            {
                elements.Add(elements[i]);
            }
            else
            {
                elements[i] = v[i - start];
            }
        }
    }

    public VectorXs normalized()
    {
        double sum = 0;
        VectorXs vectors = new VectorXs(Size);
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

    public VectorXs cross(VectorXs v)
    {
        if (v.Size != 3)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");
        VectorXs result = new VectorXs(3);
        result[0] = elements[1] * v[2] - elements[2] * v[1];
        result[1] = elements[2] * v[0] - elements[0] * v[2];
        result[2] = elements[0] * v[1] - elements[1] * v[0];
        return result;
    }

    public Vectors ToVectors()
    {
        Vectors result = new Vectors(Size);
        for (int i = 0; i < Size; i++)
        {
            result[i] = elements[i];
        }
        return result;
    }

    public VectorXs Clone()
    {
        VectorXs result = new VectorXs(Size);
        for (int i = 0; i < Size; i++)
        {
            result[i] = elements[i];
        }

        return result;
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

    public int this[int index]
    {
        get => elements[index];
        set => elements[index] = value;
    }

    public static VectorXi operator +(VectorXi a, VectorXi b)
    {
        if (a.Size != b.Size)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");

        VectorXi result = new VectorXi(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    public void SetSegment(int start, int n, VectorXi v)
    {
        if (v.Size != n)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");
        for (int i = start; i < start + n; i++)
        {
            if (i >= elements.Count)
            {
                elements.Add(elements[i]);
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
            // ũ�Ⱑ �����ϴ� ���, ������ �κ��� 0���� ä��
            elements.AddRange(new int[size - currentSize]);
        }
        else if (size < currentSize)
        {
            // ũ�Ⱑ �����ϴ� ���, ��Ҹ� �߶�
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

//��� ��Ŀ� ���
public class Triplets
{
    private int row { get; set; }
    private int col { get; set; }
    private double value { get; set; }

    public Triplets(int  r, int c, double v)
    {
        row = r;
        col = c;
        value = v;
    }
}

public class TripletXs
{
    private List<Triplets> list;

    public TripletXs()
    {
        list = new List<Triplets>();
    }

    public Triplets this[int index]
    {
        get => list[index];
        set => list[index] = value;
    }

    public void Add(int row, int col, int value)
    {
        Triplets triplets = new Triplets(row, col, value);
        list.Add(triplets);
    }

    public void Add(Triplets triplets)
    {
        list.Add(triplets);
    }
}

public class Vectors //<DIM>���� �����Ǿ� �ִ°� Ŭ���� ������ �ٲ�
{
    public int DIM;
    private double[] values;

    public Vectors(int dim)
    {
        DIM = dim;
        values = new double[dim];
    }

    public double this[int index]
    {
        get => values[index];
        set => values[index] = value;
    }

    public static Vectors operator +(Vectors a, Vectors b)
    {
        if (a.DIM != b.DIM)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");

        Vectors result = new Vectors(a.DIM);
        for (int i = 0; i < a.DIM; i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    public static Vectors operator -(Vectors a, Vectors b)
    {
        if (a.DIM != b.DIM)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");

        Vectors result = new Vectors(a.DIM);
        for (int i = 0; i < a.DIM; i++)
        {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    public static Vectors operator *(Vectors a, double b)
    {
        Vectors result = new Vectors(a.DIM);
        for (int i = 0; i < a.DIM; i++)
        {
            result[i] = a[i] * b;
        }
        return result;
    }

    public static Vectors operator *(double b, Vectors a)
    {
        Vectors result = new Vectors(a.DIM);
        for (int i = 0; i < a.DIM; i++)
        {
            result[i] = a[i] * b;
        }
        return result;
    }

    public static Vectors operator /(Vectors a, double b)
    {
        Vectors result = new Vectors(a.DIM);
        for (int i = 0; i < a.DIM; i++)
        {
            result[i] = a[i] / b;
        }
        return result;
    }

    public Vectors cross(Vectors v)
    {
        if (v.DIM != 3)
            throw new System.Exception("������ ũ�Ⱑ �ٸ��ϴ�.");
        Vectors result = new Vectors(3);
        result[0] = values[1] * v[2] - values[2] * v[1];
        result[1] = values[2] * v[0] - values[0] * v[2];
        result[2] = values[0] * v[1] - values[1] * v[0];
        return result;
    }

    public double dot(Vectors a)
    {
        double result = 0;
        for (int i = 0; i < DIM; i++)
        {
            result += values[i] * a[i];
        }
        return result;
    }

    public void setZero()
    {
        for (int i = 0; i < DIM; i++)
            values[i] = 0;
    }

    public double norm()
    {
        double result = 0;
        for (int i = 0; i < DIM; i++)
        {
            result += values[i] * values[i];
        }
        result = Math.Sqrt(result);
        return result;
    }

    public Vectors normalized()
    {
        double sum = 0;
        Vectors vectors = new Vectors(DIM);
        for (int i = 0; i < DIM; i++)
        {
            sum += values[i] * values[i];
        }

        if (sum > 0)
        {
            sum = Math.Sqrt(sum);
            for (int i = 0; i < DIM; i++)
            {
                if (values[i] != 0)
                    vectors[i] = values[i] / sum;
            }
        }
        return vectors;
    }

    public Vectors segment(int start, int n)
    {
        Vectors vectors = new Vectors(n);
        for (int i = start; i < start + n; i++)
        {
            vectors[i - start] = values[i];
        }
        return vectors;
    }

    public Vectors Clone()
    {
        Vectors result = new Vectors(DIM);
        for (int i = 0; i < DIM; i++)
        {
            result[i] = values[i];
        }
        return result;
    }

    public double squaredNorm()
    {
        double sum = 0;
        for (int i = 0; i < DIM; i++)
        {
            sum += values[i] * values[i];
        }
        return sum;
    }

    public Vectors normalize()
    {
        double sum = 0;
        Vectors vectors = new Vectors(DIM);
        for (int i = 0; i < DIM; i++)
        {
            sum += values[i] * values[i];
        }

        if (sum > 0)
        {
            sum = Math.Sqrt(sum);
            for (int i = 0; i < DIM; i++)
            {
                if (values[i] != 0)
                    vectors[i] = values[i] / sum;
            }
        }
        return vectors;
    }

    public Vectors setConstant(double value)
    {
        for (int i = 0;i < DIM;i++)
        {
            values[i] = value;
        }
        return this;
    }

    public VectorXs ToVectorXs()
    {
        VectorXs result = new VectorXs(DIM);
        for (int i = 0; i < DIM; i++)
        {
            result[i] = values[i];
        }
        return result;
    }

}

public class MatrixXs
{
    private int rows; // �� ����
    private int cols; // �� ����
    private double[,] data; // ��� ������

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
        for (int i = 0;i < rows;i++)
            for (int j = 0;j < cols;j++)
                data[i, j] = other.data[i, j];
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
            throw new IndexOutOfRangeException("�ε����� ������ ������ϴ�.");
        }
    }

    // ��Ҹ� �������� �޼���
    public double GetElement(int row, int col)
    {
        if (row >= 0 && row < rows && col >= 0 && col < cols)
        {
            return data[row, col];
        }
        else
        {
            throw new IndexOutOfRangeException("�ε����� ������ ������ϴ�.");
        }
    }

    public void setZero()
    {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                data[i, j] = 0;
    }

    // ��� ���� �޼���
    public static MatrixXs operator+(MatrixXs a, MatrixXs b)
    {
        if (a.rows != b.rows || a.cols != b.cols)
        {
            throw new InvalidOperationException("����� ũ�Ⱑ ��ġ���� �ʽ��ϴ�.");
        }

        MatrixXs result = new MatrixXs(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
            {
                result.SetElement(i, j, a.GetElement(i, j) + b.GetElement(i, j));
            }
        }

        return result;
    }

    // ��� ���� �޼���
    public static MatrixXs operator*(MatrixXs a, MatrixXs b)
    {
        if (a.cols != b.rows)
        {
            throw new InvalidOperationException("ù ��° ����� ���� ������ �� ��° ����� ���� ������ ��ġ�ؾ� �մϴ�.");
        }

        MatrixXs result = new MatrixXs(a.rows, b.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < b.cols; j++)
            {
                double sum = 0;
                for (int k = 0; k < a.cols; k++)
                {
                    sum += a.GetElement(i, k) * b.GetElement(k, j);
                }
                result.SetElement(i, j, sum);
            }
        }
        return result;
    }

    public static MatrixXs operator*(double a, MatrixXs b)
    {
        MatrixXs result = new MatrixXs(b.rows, b.cols);
        for (int i = 0; i < b.rows; i++)
        {
            for (int j = 0; j < b.cols; j++)
            {
                result[i, j] = b[i, j];
                result[i, j] *= a;
            }
        }
        return result;
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
                matrix.data[j, i] = data[i, j];
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

    public VectorXs col(int i)
    {
        VectorXs v = new VectorXs(i);
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
    private double angle;  // ȸ�� ���� (����)

    public Rotation2D(double angleInRadians)
    {
        angle = angleInRadians;
    }

    // ȸ�� ����� ����Ͽ� ���� (x, y)�� ȸ��
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
    private Dictionary<(int, int), double> matrix;
    public int Rows { get; private set; }
    public int Cols { get; private set; }

    public SparseXs(int rows, int cols)
    {
        matrix = new Dictionary<(int, int), double>();
        Rows = rows;
        Cols = cols;
    }

    public void insert(int row, int col, double value)
    {
        if (row >= Rows || col >= Cols || row < 0 || col < 0)
        {
            throw new ArgumentOutOfRangeException("Invalid matrix indices.");
        }

        if (value == 0)
        {
            matrix.Remove((row, col));
        }
        else
        {
            matrix[(row, col)] = value;
        }
    }

    public double GetValue(int row, int col)
    {
        if (row >= Rows || col >= Cols || row < 0 || col < 0)
        {
            throw new ArgumentOutOfRangeException("Invalid matrix indices.");
        }

        return matrix.ContainsKey((row, col)) ? matrix[(row, col)] : 0.0;
    }

}