using UnityEngine;
using UnityEngine.Rendering;

[RequireComponent(typeof(Camera))]
public class MRTScript : MonoBehaviour
{
    public RenderTexture colorRenderTexture;   // 색상 정보 저장용
    public RenderTexture normalRenderTexture;  // 노멀 정보 저장용
    public RenderTexture reflectUvRenderTexture; // SSR 반사 정보 저장용
    public Renderer target;

    public Camera cam; // 물을 렌더링할 카메라. 메인 카메라면 안 됨
    private Material mrtMaterial;
    private CommandBuffer commandBuffer;

    private void Start()
    {
        // MRT용 Material 설정
        mrtMaterial = target.material;
        
        // 커맨드 버퍼 초기화
        commandBuffer = new CommandBuffer();
        commandBuffer.name = "MRT Command Buffer";

        // MRT 설정 - 여러 개의 Render Target에 출력 설정
        RenderTargetIdentifier[] mrt = new RenderTargetIdentifier[3]
        {
            colorRenderTexture,
            normalRenderTexture,
            reflectUvRenderTexture,
        };

        // 두 Render Texture에 MRT 셰이더로 렌더링
        commandBuffer.SetRenderTarget(mrt, colorRenderTexture.depthBuffer);
        commandBuffer.ClearRenderTarget(true, true, Color.clear);
        commandBuffer.DrawRenderer(target, mrtMaterial);

        // 카메라 렌더링 직전에 커맨드 버퍼 추가
        cam.AddCommandBuffer(CameraEvent.BeforeForwardOpaque, commandBuffer);
        //cam.AddCommandBuffer(CameraEvent.AfterEverything, commandBuffer);
    }

    private void OnDisable()
    {
        // 커맨드 버퍼 해제
        if (cam != null)
        {
            //cam.RemoveCommandBuffer(CameraEvent.BeforeForwardOpaque, commandBuffer);
        }
    }
}
