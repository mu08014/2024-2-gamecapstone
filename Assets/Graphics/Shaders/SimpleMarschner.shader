Shader "Custom/SimpleMarschner"
{
    Properties
    {
        _DiffuseTerm("DiffuseTerm", Range(0.0, 1.0)) = 1
        _SpecularTerm("SpecularTerm", Range(0.0, 1.0)) = 1
        _SpecularPower("SpecularPower", Range(0.0, 1000.0)) = 10
        _Color ("Color", Color) = (1,1,1,1)
        _SpecularColor0 ("SpecularColor0", Color) = (1,1,1,1)
        _SpecularColor1 ("SpecularColor1", Color) = (1,1,1,1)
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        //LOD 200
        
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #include "UnityCG.cginc"
            
            float _DiffuseTerm;
            float _SpecularTerm;
            float _SpecularPower;
            fixed4 _Color;
            float4 _SpecularColor0;
            float4 _SpecularColor1;
            sampler2D _MainTex;
            
            struct appdata
            {
                float4 vertex : POSITION;
                float3 tangent : TEXCOORD1;
            };
            
            struct v2f
            {
                //float2 uv : TEXCOORD0;
                float4 pos : SV_POSITION;
                float3 tangent : TANGENT;
                float3 lightDir : TEXCOORD1;
                float3 cameraDir : TEXCOORD2;
            };
            
            v2f vert (appdata v)
            {
                v2f o;
                o.pos = UnityObjectToClipPos(v.vertex);
                o.tangent = normalize(v.tangent);
                o.lightDir = mul(unity_WorldToObject, _WorldSpaceLightPos0.xyz);
                o.cameraDir = normalize(
                    mul(unity_WorldToObject, float4(_WorldSpaceCameraPos.xyz, 1.0f))
                    - v.vertex); 
                return o;
            }
            
            float4 frag (v2f i) : SV_Target
            {
                float3 diffuse = _Color * _DiffuseTerm;

                float3 halfVec = normalize(i.lightDir + i.cameraDir);
                
                float3 shift_tangent0 = normalize(i.tangent + float3(0, 1, 0) * 0.25f);
                float3 shift_tangent1 = normalize(i.tangent - float3(0, 1, 0) * 0.25f);

                float3 specular0 = pow( 
                    sqrt(1 - dot(shift_tangent0, halfVec) * dot(shift_tangent0, halfVec)),
                    _SpecularPower) * _SpecularTerm;
                
                float3 specular1 = pow( 
                    sqrt(1 - dot(shift_tangent1, halfVec) * dot(shift_tangent1, halfVec)),
                    _SpecularPower) * _SpecularTerm;
                  

                float3 color = diffuse + specular1 * _SpecularColor0 + specular1 * _SpecularColor1;
                
                return float4(color, 1);
            }
            ENDCG
        }
    }
    FallBack "Diffuse"
}
