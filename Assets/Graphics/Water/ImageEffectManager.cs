using Unity.VisualScripting;
using UnityEngine;

[ExecuteInEditMode]
public class ImageEffectManager : MonoBehaviour
{
    public Material imageEffectMaterial; // ���̴��� ����� Material
    public Texture distortionTex = null; // �ְ� �ؽ�ó

    private void OnRenderImage(RenderTexture src, RenderTexture dest)
    {
        if (imageEffectMaterial != null && distortionTex != null)
        {
            // ���̴��� ����� Material�� ����� src�� dest�� ������
            imageEffectMaterial.SetTexture("_NoiseTex", distortionTex);
            Graphics.Blit(src, dest, imageEffectMaterial);
        }
    }
}
