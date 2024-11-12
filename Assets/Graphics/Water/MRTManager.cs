using UnityEngine;
using UnityEngine.Rendering;

[RequireComponent(typeof(Camera))]
public class MRTScript : MonoBehaviour
{
    public RenderTexture colorRenderTexture;   // ���� ���� �����
    public RenderTexture normalRenderTexture;  // ��� ���� �����
    public Renderer target;

    private Camera cam;
    private Material mrtMaterial;
    private CommandBuffer commandBuffer;

    private void Start()
    {
        //colorRenderTexture = new(1920, 1080, 1);
        //normalRenderTexture = new(1920, 1080, 1);

        cam = GetComponent<Camera>();

        
        // MRT�� Material ����
        mrtMaterial = target.material;
        
        // Ŀ�ǵ� ���� �ʱ�ȭ
        commandBuffer = new CommandBuffer();
        commandBuffer.name = "MRT Command Buffer";

        // MRT ���� - ���� ���� Render Target�� ��� ����
        RenderTargetIdentifier[] mrt = new RenderTargetIdentifier[2]
        {
            //normalRenderTexture,
            colorRenderTexture,
            normalRenderTexture
        };

        // �� Render Texture�� MRT ���̴��� ������
        commandBuffer.SetRenderTarget(mrt, colorRenderTexture.depthBuffer);
        commandBuffer.ClearRenderTarget(true, true, Color.clear);
        commandBuffer.DrawRenderer(target, mrtMaterial);

        // ī�޶� ������ ������ Ŀ�ǵ� ���� �߰�
        cam.AddCommandBuffer(CameraEvent.BeforeForwardOpaque, commandBuffer);
    }

    private void OnDisable()
    {
        // Ŀ�ǵ� ���� ����
        if (cam != null)
        {
            cam.RemoveCommandBuffer(CameraEvent.BeforeForwardOpaque, commandBuffer);
        }
    }
}
