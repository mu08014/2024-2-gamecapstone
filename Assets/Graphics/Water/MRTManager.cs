using UnityEngine;
using UnityEngine.Rendering;

[RequireComponent(typeof(Camera))]
public class MRTScript : MonoBehaviour
{
    public RenderTexture colorRenderTexture;   // ���� ���� �����
    public RenderTexture normalRenderTexture;  // ��� ���� �����
    public RenderTexture reflectUvRenderTexture; // SSR �ݻ� ���� �����
    public Renderer target;

    public Camera cam; // ���� �������� ī�޶�. ���� ī�޶�� �� ��
    private Material mrtMaterial;
    private CommandBuffer commandBuffer;

    private void Start()
    {
        // MRT�� Material ����
        mrtMaterial = target.material;
        
        // Ŀ�ǵ� ���� �ʱ�ȭ
        commandBuffer = new CommandBuffer();
        commandBuffer.name = "MRT Command Buffer";

        // MRT ���� - ���� ���� Render Target�� ��� ����
        RenderTargetIdentifier[] mrt = new RenderTargetIdentifier[3]
        {
            colorRenderTexture,
            normalRenderTexture,
            reflectUvRenderTexture,
        };

        // �� Render Texture�� MRT ���̴��� ������
        commandBuffer.SetRenderTarget(mrt, colorRenderTexture.depthBuffer);
        commandBuffer.ClearRenderTarget(true, true, Color.clear);
        commandBuffer.DrawRenderer(target, mrtMaterial);

        // ī�޶� ������ ������ Ŀ�ǵ� ���� �߰�
        cam.AddCommandBuffer(CameraEvent.BeforeForwardOpaque, commandBuffer);
        //cam.AddCommandBuffer(CameraEvent.AfterEverything, commandBuffer);
    }

    private void OnDisable()
    {
        // Ŀ�ǵ� ���� ����
        if (cam != null)
        {
            //cam.RemoveCommandBuffer(CameraEvent.BeforeForwardOpaque, commandBuffer);
        }
    }
}
