using UnityEngine;
using UnityEngine.Rendering;

[RequireComponent(typeof(Camera))]
public class MRTScript : MonoBehaviour
{
    public RenderTexture colorRenderTexture;   // 색상 정보 저장용
    public RenderTexture normalRenderTexture;  // 노멀 정보 저장용
    public Renderer target;

    private Camera cam;
    private Material mrtMaterial;
    private CommandBuffer commandBuffer;

    private void Start()
    {
        //colorRenderTexture = new(1920, 1080, 1);
        //normalRenderTexture = new(1920, 1080, 1);

        cam = GetComponent<Camera>();

        
        // MRT용 Material 설정
        mrtMaterial = target.material;
        
        // 커맨드 버퍼 초기화
        commandBuffer = new CommandBuffer();
        commandBuffer.name = "MRT Command Buffer";

        // MRT 설정 - 여러 개의 Render Target에 출력 설정
        RenderTargetIdentifier[] mrt = new RenderTargetIdentifier[2]
        {
            //normalRenderTexture,
            colorRenderTexture,
            normalRenderTexture
        };

        // 두 Render Texture에 MRT 셰이더로 렌더링
        commandBuffer.SetRenderTarget(mrt, colorRenderTexture.depthBuffer);
        commandBuffer.ClearRenderTarget(true, true, Color.clear);
        commandBuffer.DrawRenderer(target, mrtMaterial);

        // 카메라 렌더링 직전에 커맨드 버퍼 추가
        cam.AddCommandBuffer(CameraEvent.BeforeForwardOpaque, commandBuffer);
    }

    private void OnDisable()
    {
        // 커맨드 버퍼 해제
        if (cam != null)
        {
            cam.RemoveCommandBuffer(CameraEvent.BeforeForwardOpaque, commandBuffer);
        }
    }
}
