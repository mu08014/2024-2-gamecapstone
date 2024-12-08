using Unity.VisualScripting;
using UnityEngine;

// Add water. 
public class ImageEffectManager : MonoBehaviour
{
    public Material imageEffectMaterial; // 셰이더가 적용된 Material
    public Texture colorTex = null; // 색 텍스처
    public Texture distortionTex = null; // 왜곡 텍스처

    private void Start()
    {
        imageEffectMaterial.SetTexture("_ColorTex", colorTex);
        imageEffectMaterial.SetTexture("_NoiseTex", distortionTex);
    }

    private void OnRenderImage(RenderTexture src, RenderTexture dest)
    {
        if (imageEffectMaterial != null && distortionTex != null)
        {
            // 셰이더가 적용된 Material을 사용해 src를 dest로 렌더링
            Graphics.Blit(src, dest, imageEffectMaterial);
        }
    }
}
