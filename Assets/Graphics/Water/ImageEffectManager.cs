using UnityEngine;

[ExecuteInEditMode]
public class ImageEffectManager : MonoBehaviour
{
    public Material imageEffectMaterial; // 셰이더가 적용된 Material

    private void OnRenderImage(RenderTexture src, RenderTexture dest)
    {
        if (imageEffectMaterial != null)
        {
            // 셰이더가 적용된 Material을 사용해 src를 dest로 렌더링
            Graphics.Blit(src, dest, imageEffectMaterial);
        }
    }
}
