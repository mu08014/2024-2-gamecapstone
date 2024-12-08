using Unity.VisualScripting;
using UnityEngine;

// Add water. 
public class ImageEffectManager : MonoBehaviour
{
    public Material imageEffectMaterial; // ���̴��� ����� Material
    public Texture colorTex = null; // �� �ؽ�ó
    public Texture distortionTex = null; // �ְ� �ؽ�ó

    private void Start()
    {
        imageEffectMaterial.SetTexture("_ColorTex", colorTex);
        imageEffectMaterial.SetTexture("_NoiseTex", distortionTex);
    }

    private void OnRenderImage(RenderTexture src, RenderTexture dest)
    {
        if (imageEffectMaterial != null && distortionTex != null)
        {
            // ���̴��� ����� Material�� ����� src�� dest�� ������
            Graphics.Blit(src, dest, imageEffectMaterial);
        }
    }
}
