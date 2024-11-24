using UnityEngine;

public class GLTest : MonoBehaviour
{
    // When added to an object, draws colored rays from the
    // transform position.
    public int lineCount = 100;
    public float length = 1.0f;

    //https://docs.unity3d.com/ScriptReference/GL.html
    public Material lineMaterial;
    private void CreateLineMaterial()
    {
        if (!lineMaterial)
        {
            // Unity has a built-in shader that is useful for drawing
            // simple colored things.
            Shader shader = Shader.Find("Hidden/Internal-Colored");
            lineMaterial = new Material(shader);
            lineMaterial.hideFlags = HideFlags.HideAndDontSave;
            // Turn on alpha blending
            lineMaterial.SetInt("_SrcBlend", (int)UnityEngine.Rendering.BlendMode.SrcAlpha);
            lineMaterial.SetInt("_DstBlend", (int)UnityEngine.Rendering.BlendMode.OneMinusSrcAlpha);
            // Turn backface culling off
            lineMaterial.SetInt("_Cull", (int)UnityEngine.Rendering.CullMode.Off);
            // Turn off depth writes
            lineMaterial.SetInt("_ZWrite", 0);
        }
    }

    // Will be called after all regular rendering is done
    public void OnRenderObject()
    {
        CreateLineMaterial();
        // Apply the line material
        lineMaterial.SetPass(0);

        GL.PushMatrix();
        // Set transformation matrix for drawing to
        // match our transform
        GL.MultMatrix(transform.localToWorldMatrix);

        // Draw lines
        GL.Begin(GL.LINE_STRIP);
        for (int i = 0; i < lineCount; ++i)
        {
            for (int j = 0; j < lineCount; ++j)
            {
                Vector3 pos = new Vector3(1.0f / lineCount * i, 0.0f, 1.0f / lineCount * j);

                GL.Color(Color.red);

                GL.Vertex3(pos.x, pos.y, pos.z);
                GL.Vertex3(pos.x, pos.y + length, pos.z);
            }    
        }

        GL.End();
        GL.PopMatrix();
    }
}
