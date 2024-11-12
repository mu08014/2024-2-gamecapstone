using UnityEngine;

[ExecuteInEditMode]
public class ImageEffectManager : MonoBehaviour
{
    public Material imageEffectMaterial; // ���̴��� ����� Material

    private void OnRenderImage(RenderTexture src, RenderTexture dest)
    {
        if (imageEffectMaterial != null)
        {
            // ���̴��� ����� Material�� ����� src�� dest�� ������
            Graphics.Blit(src, dest, imageEffectMaterial);
        }
    }
}
