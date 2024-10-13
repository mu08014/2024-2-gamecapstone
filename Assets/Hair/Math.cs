using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using System.Xml.Linq;
using Unity.VisualScripting;
using UnityEngine;

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
}

public class VectorXs
{
    private List<double> elements;
    private int Size => elements.Count;

    public VectorXs()
    {
        elements = new List<double>();
    }

    public VectorXs(int size)
    {
        elements = new List<double>(size);
        for (int i = 0; i < size; i++)
        {
            elements[i] = 0;
        }
    }

    public double this[int index]
    {
        get => elements[index];
        set => elements[index] = value;
    }

    public static VectorXs operator +(VectorXs a, VectorXs b)
    {
        if (a.Size != b.Size)
            throw new System.Exception("벡터의 크기가 다릅니다.");

        VectorXs result = new VectorXs(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    public int size()
    {
        return elements.Count;
    }

    public void resize(int size)
    {
        List<double> values = new List<double>(size);
        for (int i = 0; i < size; i++)
        {
            values[i] = elements[i];
        }

        elements = new List<double>(size);
        for (int i = 0; i < size; i++)
        {
            elements[i] = values[i];
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

    // .segment() = Vectors 형태일 때 사용
    public void SetSegment(int start, int n, Vectors v)
    {
        if (v.DIM != n)
            throw new System.Exception("벡터의 크기가 다릅니다.");

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
            throw new System.Exception("벡터의 크기가 다릅니다.");

        VectorXi result = new VectorXi(a.Size);
        for (int i = 0; i < a.Size; i++)
        {
            result[i] = a[i] + b[i];
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

public class Vectors //<DIM>으로 구현되어 있는거 클래스 변수로 바꿈
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
            throw new System.Exception("벡터의 크기가 다릅니다.");

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
            throw new System.Exception("벡터의 크기가 다릅니다.");

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
                vectors[i] = values[i] * sum;
            }
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

    public void setZero()
    {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                data[i, j] = 0;
    }

    // 행렬 덧셈 메서드
    public static MatrixXs operator+(MatrixXs a, MatrixXs b)
    {
        if (a.rows != b.rows || a.cols != b.cols)
        {
            throw new InvalidOperationException("행렬의 크기가 일치하지 않습니다.");
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

    // 행렬 곱셈 메서드
    public static MatrixXs operator*(MatrixXs a, MatrixXs b)
    {
        if (a.cols != b.rows)
        {
            throw new InvalidOperationException("첫 번째 행렬의 열의 개수와 두 번째 행렬의 행의 개수가 일치해야 합니다.");
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