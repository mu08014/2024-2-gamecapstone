using Unity.VisualScripting;
using UnityEngine;

[ExecuteInEditMode]
public class ImageEffectManager : MonoBehaviour
{
    public Material imageEffectMaterial; // 셰이더가 적용된 Material
    public Texture distortionTex = null; // 왜곡 텍스처

    private void OnRenderImage(RenderTexture src, RenderTexture dest)
    {
        if (imageEffectMaterial != null && distortionTex != null)
        {
            // 셰이더가 적용된 Material을 사용해 src를 dest로 렌더링
            imageEffectMaterial.SetTexture("_NoiseTex", distortionTex);
            Graphics.Blit(src, dest, imageEffectMaterial);
        }
    }
}
